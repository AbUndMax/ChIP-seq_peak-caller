import itertools
import threading
import time
import math
import sys


# self-made module for loading-bars and throbbers


class Throbber:
    def __init__(self, desc="Loading...", end="Done!", timeout=0.3):
        self.desc = desc
        self.end = end
        self.timeout = timeout
        self._running = False
        self._thread = None

    def start(self):
        self._running = True
        self._thread = threading.Thread(target=self._animate)
        self._thread.start()

    def _animate(self):
        for c in itertools.cycle(['|', '/', '-', '\\']):
            if not self._running:
                break
            sys.stdout.write(f'\r{self.desc} {c}')
            sys.stdout.flush()
            time.sleep(self.timeout)
        sys.stdout.write(f'\r{self.end}\n')
        sys.stdout.flush()

    def stop(self):
        self._running = False
        if self._thread is not None:
            self._thread.join()


class LoadingBar:
    bar_end = "] % finished"

    def __init__(self, total, desc="Loading"):
        self.total = total - 1
        self.finished = False
        self.bar_start = desc + ": ["
        self.__percent = -1
        self.load(0)
        self.__percent = 0

    def load(self, current):
        current_percent = int(current / self.total * 100)

        if self.__percent < current_percent < 100:
            self.__percent = current_percent
            self._print_bar(current_percent)

        elif not self.finished and current_percent == 100:
            self._print_bar(100, final=True)

        else:
            return

    def _print_bar(self, current_percent, final=False):
        current_percent_string = f" {current_percent:02.0f}"
        repetitions = int(math.floor(math.floor(current_percent) / 10.0))

        bar = self.bar_start
        bar = bar + " ##" * (repetitions if repetitions < 9 else 9) + current_percent_string
        bar += " --" * (10 - repetitions - 1) + ("" if current_percent >= 100 else " ") + self.bar_end

        if final:
            self.finished = True
            sys.stdout.write(f'\r{bar}\n\n')
        else:
            sys.stdout.write(f'\r{bar}')

        sys.stdout.flush()
