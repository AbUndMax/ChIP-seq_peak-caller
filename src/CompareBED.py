import matplotlib.pyplot as plt
from matplotlib_venn import venn2


class Peak:
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end

    def overlaps(self, other):
        return self.chrom == other.chrom and self.start < other.end and self.end > other.start


def read_bed(file):
    peaks = []
    with open(file, 'r') as f:
        for line in f:
            chrom, start, end = line.strip().split('\t')[:3]
            peaks.append(Peak(chrom, int(start), int(end)))
    return peaks


def find_overlaps(peaks1, peaks2):
    overlaps = set()
    for p1 in peaks1:
        for p2 in peaks2:
            if p1.overlaps(p2):
                overlaps.add(p1)
                overlaps.add(p2)
    return overlaps

# Datei-Pfade
MACS_peaks = '/Users/max/Library/CloudStorage/OneDrive-Persönlich/Bioinformatics_M.Sc/Module/2_Semester_M.Sc/Next_Generation_Sequencing_with_Python/NGS_assignments/Assignment_07_Files/ChIP-seq/Oct4_summits.bed'
MY_peaks = '/Users/max/Library/CloudStorage/OneDrive-Persönlich/Bioinformatics_M.Sc/Module/2_Semester_M.Sc/Next_Generation_Sequencing_with_Python/NGS_assignments/Assignment_07_Files/asgmt07_NiklasGerbes_6685993/called_peaks.bed'

# Peaks einlesen
my_peaks = read_bed(MY_peaks)
macs_peaks = read_bed(MACS_peaks)

# Überlappende Peaks finden
overlapping_peaks = find_overlaps(my_peaks, macs_peaks)

# Mengen für Venn-Diagramm vorbereiten
my_peaks_set = set(my_peaks)
macs_peaks_set = set(macs_peaks)
overlapping_peaks_set = set(overlapping_peaks)

# Größen für das Venn-Diagramm berechnen
only_my_peaks = len(my_peaks_set - overlapping_peaks_set)
only_macs_peaks = len(macs_peaks_set - overlapping_peaks_set)
common_peaks = len(overlapping_peaks_set)

# Venn-Diagramm plotten
plt.figure(figsize=(8, 8))
venn2(subsets=(only_my_peaks, only_macs_peaks, common_peaks), set_labels=('My Peaks', 'MACS Peaks'))
plt.title('Overlap of Peaks from Two BED Files')
plt.show()