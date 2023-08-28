import time


class time_meter:
    def __init__(self):
        self.s = time.time_ns()

    def meter(self):
        t = time.time_ns()
        return (t - self.s) / 1e6
