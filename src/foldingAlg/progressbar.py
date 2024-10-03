import sys
import time

# simple progress bar class, mostly chatgpt code, 2024

class ProgressBar:
    """
    Progress bar initialization example
        progress_bar = ProgressBar(total=10, prefix="Processing: ")
        for _ in range(10):
            time.sleep(0.5)
            progress_bar.update()

    initialization with iterator
        dummy_function = lambda: [time.sleep(0.5) for _ in ProgressBar.iterate(range(4), "Processing: ")]
        dummy_function()   
    """
        
    def __init__(self, total, prefix="", size=60, out=sys.stdout):
        self.total = total
        self.current = 0
        self.prefix = prefix
        self.size = size
        self.out = out
        self.start = time.time()
        self.full_block = "█"
        self.three_quarters_block = "▊"
        self.half_block = "▌"
        self.quarter_block = "▎"
        self.no_block = "."

    def update(self, increment=1):
        self.current = min(self.current + increment, self.total)
        self._display_bar()

    def _display_bar(self):
        progress = self.current / self.total
        x_full = int(self.size * progress)
        bar = self.full_block * x_full + self.no_block * (self.size - x_full)

        if self.current < self.total:
            remaining = ((time.time() - self.start) / max(1, self.current)) * (self.total - self.current)
            mins, sec = divmod(remaining, 60)
            time_str = f"{int(mins):02d}:{sec:05.2f}"
            print(f"{self.prefix}[{bar}] {self.current}/{self.total} Est wait {time_str}", end='\r', file=self.out, flush=True)
        else:
            self.complete()

    def reset(self):
        self.current = 0
        self.start = time.time()


    def complete(self):
        total_time = time.time() - self.start
        mins, sec = divmod(total_time, 60)
        # time_str = f"{int(mins)} min {int(sec)} sec"
        time_str = f"{int(mins):02d}:{sec:05.2f}"
        bar = self.full_block * self.size
        print(f"{self.prefix}[{bar}] {self.total}/{self.total} Completed in {time_str}", end='\n', file=self.out, flush=True)


    @classmethod
    def iterate(cls, iterable, prefix="", size=60, out=sys.stdout):
        progress_bar = cls(total=len(iterable), prefix=prefix, size=size, out=out)
        for item in iterable:
            yield item
            progress_bar.update()