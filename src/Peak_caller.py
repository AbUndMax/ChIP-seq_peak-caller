import pysam
from scipy.stats import poisson
from numpy import mean, log10
from pathlib import Path
import argparse
import loading
import os


# this class is used as a data structure to store the information of a peak
class Peak:

    def __init__(self, chromosome, direction, start, end, peak_coverage, peak_max_position):
        self.chromosome = chromosome
        self.direction = direction
        self.start = start
        self.end = end
        self.coverage = peak_coverage
        self.peak_max_position = peak_max_position
        self.length = end - start + 1
        self.p_value = 0.0

    def __str__(self):
        return f"Peak start: {self.start}, Peak end: {self.end}, Peak coverage: {self.coverage}, Peak length: {self.length}, Peak max position: {self.peak_max_position}"

    def to_BED_row(self):
        return f"{self.chromosome}\t{self.start}\t{self.end}\t{self.p_value}"


# Checks if the index file (.bai) for the BAM file exists
def check_bam_index(bam_path):
    bam_file = Path(bam_path)
    bai_file = bam_file.with_suffix(bam_file.suffix + '.bai')
    
    if bai_file.exists():
        return True
    else:
        print(f"Missing .bai file: {bam_file.name}")
        print("> Make sure the .bai is named the same as the bam file and is in the same directory!")
        return False
    

# Calculate coverage per chromosome, direction and position -> output is dictioniary with chromsome as key and a new
# dictionary as value. This new dictionary has the direction as key and another dictionary as value. This last dictionary
# has the position as key and the coverage as value.
# the encapsulation of this nested dictionary helps to loop over the chromosomes and directions
def calculate_coverage(bam_file):
    throbbing = loading.Throbber(desc="Calculating coverage")
    throbbing.start()

    bam = pysam.AlignmentFile(bam_file, "rb")

    coverage_per_chromosome_dict = {}

    for chromosome in bam.references:  # Iterate through each chromosome
        coverage_per_chromosome_dict[chromosome] = {
            'forward': {},
            'reverse': {}
        }

        for pileup_column in bam.pileup(contig=chromosome):  # Specify the chromosome for the pileup
            pos = pileup_column.reference_pos
            for pileup_read in pileup_column.pileups:
                if not pileup_read.is_del and not pileup_read.is_refskip:
                    # if the read is reverse, add the coverage to the reverse dictionary, else to the forward dictionary
                    # then increment the read count by one
                    if pileup_read.alignment.is_reverse:
                        if pos in coverage_per_chromosome_dict[chromosome]['reverse']:
                            coverage_per_chromosome_dict[chromosome]['reverse'][pos] += 1
                        else:
                            coverage_per_chromosome_dict[chromosome]['reverse'][pos] = 1
                    else:
                        if pos in coverage_per_chromosome_dict[chromosome]['forward']:
                            coverage_per_chromosome_dict[chromosome]['forward'][pos] += 1
                        else:
                            coverage_per_chromosome_dict[chromosome]['forward'][pos] = 1

    bam.close()
    throbbing.stop()
    return coverage_per_chromosome_dict


# Calculate the average coverage of the background (average coverage equals the expected value in poisson)
def calculate_average_background_coverage(coverage_per_chromosome_dict):
    all_coverages = []
    for chromosome in coverage_per_chromosome_dict:
        all_coverages += list(coverage_per_chromosome_dict[chromosome]['forward'].values())
        all_coverages += list(coverage_per_chromosome_dict[chromosome]['reverse'].values())

    average_coverage = mean(all_coverages)
    print(f"Average coverage: {average_coverage}\n")
    return average_coverage


# I use an initial window to search for high coverage regions. If a window has a mean coverage above the average
# coverage, I extend the window to include adjacent positions with coverage above average.
def new_peak_identifier(coverage_per_chromosome_dict, average_coverage, window_size=10):
    peaks_per_chromosome_dict = {}

    # Iterate over all chromosomes
    for chromosome in coverage_per_chromosome_dict:

        peaks_per_chromosome_dict[chromosome] = {
            'forward': [],
            'reverse': []
        }

        # Iterate over all directions (reverse and forward)
        for direction in coverage_per_chromosome_dict[chromosome]:
            coverage_dict = coverage_per_chromosome_dict[chromosome][direction]
            positions = list(coverage_dict.keys())
            positions.sort()

            bar = loading.LoadingBar(len(positions) - 1,
                                     f"searching peaks on chromosome {chromosome} in {direction} direction")

            window_start = 0
            while window_start < len(positions):  # Iterate over all positions
                window_end = min(window_start + window_size, len(positions))
                window = positions[window_start:window_end]

                if window_is_continuous(window):  # Check if the positions in the window are continuous
                    window_coverage_mean = mean([coverage_dict[position] for position in window])

                    if window_coverage_mean > average_coverage:
                        # Extend the window to include adjacent positions with coverage above average

                        next_position = positions[window_end]
                        while next_position in coverage_dict and next_position == positions[window_end - 1] + 1 and \
                                coverage_dict[next_position] > average_coverage and window_end < len(positions) - 1:
                            window_end += 1
                            next_position = positions[window_end]
                            bar.load(window_end)

                        peak_coverage_per_position = [coverage_dict[pos] for pos in positions[window_start:window_end]]
                        peak_coverage = sum(peak_coverage_per_position)
                        peak_max_coverage = max(peak_coverage_per_position)

                        # peak_max_position is the position where the coverage is the last of the most highest
                        # i.e. if the read is forward directed the position of the max is the most right position
                        # if the read is reverse directed the position of the max is the most left position
                        if direction == 'forward':
                            peak_max_position = positions[
                                window_start + rindex(peak_coverage_per_position, peak_max_coverage)]
                        else:
                            peak_max_position = positions[
                                window_start + peak_coverage_per_position.index(peak_max_coverage)]
                        peaks_per_chromosome_dict[chromosome][direction].append(
                            Peak(chromosome, direction, positions[window_start], positions[window_end - 1], peak_coverage,
                                 peak_max_position))  # window_end -1 because end is exclusive
                    window_start = window_end
                else:
                    window_start += 1
                    bar.load(window_end)
                bar.load(window_end)

    return peaks_per_chromosome_dict


def rindex(lst, value):
    return len(lst) - lst[::-1].index(value) - 1


def window_is_continuous(window: list):
    return window[-1] - window[0] == len(window) - 1


def filter_all_peaks(peaks_per_chromosome, average_coverage, p_value_threshold=0.01):
    filtered_peaks_per_chromosome = {}

    # again we loop over all chromosomes and directions
    for chromosome in peaks_per_chromosome:
        filtered_peaks_per_chromosome[chromosome] = {
            'forward': [],
            'reverse': []
        }

        # filter peaks with a p-value below the threshold
        for direction in peaks_per_chromosome[chromosome]:
            peaks = peaks_per_chromosome[chromosome][direction]
            filtered_peaks = filter_peaks(peaks, average_coverage, chromosome, direction, p_value_threshold)
            filtered_peaks_per_chromosome[chromosome][direction] = filtered_peaks

    return filtered_peaks_per_chromosome


def filter_peaks(peaks, average_coverage, chromosome, direction, p_value_threshold=0.01):
    filtered_peaks = []

    bar = loading.LoadingBar(len(peaks) - 1, desc=f"Filtering peaks on {chromosome} in {direction} direction")
    for peak in peaks:
        peak_p_value = poisson.sf(peak.coverage - 1, average_coverage * peak.length)
        if peak_p_value < p_value_threshold:
            peak.p_value = -log10(peak_p_value) if peak_p_value != 0 else 0
            filtered_peaks.append(peak)
        bar.load(peaks.index(peak))

    return filtered_peaks


def call_peaks(peaks_per_chromosome, average_coverage, p_value_threshold=0.01, peak_distance=200):
    filtered_peaks = filter_all_peaks(peaks_per_chromosome, average_coverage, p_value_threshold)

    called_peaks = []
    for chromosome in filtered_peaks:
        for peak_fwd in filtered_peaks[chromosome]['forward']:
            for peak_rev in filtered_peaks[chromosome]['reverse']:
                # we search for peaks that have their max position within peak_distance of each other
                if peak_rev.peak_max_position - peak_distance < peak_fwd.peak_max_position < peak_rev.peak_max_position:
                    peak_max_delta = peak_rev.peak_max_position - peak_fwd.peak_max_position
                    # then we shift the start and end of the forward peak by half of the peak_distance
                    shifted_start = int(peak_fwd.start + peak_max_delta // 2)
                    shifted_end = int(peak_fwd.end + peak_max_delta // 2)
                    shifted_max = int(peak_fwd.peak_max_position + peak_max_delta // 2)
                    shifted_peak = Peak(peak_fwd.chromosome, peak_fwd.direction, shifted_start, shifted_end,
                                        peak_fwd.coverage, shifted_max)
                    # I just used the mean of the p-values of the two peaks
                    # I know this is not a good way to combine the p-values!
                    shifted_peak.p_value = mean([peak_fwd.p_value, peak_rev.p_value])
                    called_peaks.append(shifted_peak)
    return called_peaks


# simple function to write the peaks to a BED file
def write_BED_file(peaks, output_path):
    with open(output_path, 'w') as f:
        for peak in peaks:
            f.write(peak.to_BED_row() + '\n')


def main():
    output_file = os.getcwd() + '/called_peaks.bed'
    parser = argparse.ArgumentParser(description='Peak caller')
    parser.add_argument('--BAM_path', '-b', type=str, help='Path to the BAM file', required=True)
    parser.add_argument('--output_path', '-o', type=str, help='Path to the output file', default=output_file)
    args = parser.parse_args()
    
    bam_path = Path(args.BAM_path)
     
    if not check_bam_index(bam_path):
        return

    coverage_per_chromosome_dict = calculate_coverage(bam_path)
    average_coverage = calculate_average_background_coverage(coverage_per_chromosome_dict)
    peaks_per_chromosome = new_peak_identifier(coverage_per_chromosome_dict, average_coverage, 40)
    called_peaks = call_peaks(peaks_per_chromosome, average_coverage)
    write_BED_file(called_peaks, args.output_path)


if __name__ == '__main__':
    main()
