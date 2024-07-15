from abc import ABC, abstractmethod
from dataclasses import dataclass


class MappingFlagBuild(ABC):
    build_name = None

    def __init__(self, depth, depthc, coverage, mapped, windows_covered):
        self.depth = depth
        self.depthc = depthc
        self.coverage = coverage
        self.mapped = mapped
        self.windows_covered = windows_covered

    @abstractmethod
    def assert_false_positive(self) -> bool:
        pass

    @abstractmethod
    def assert_vestigial(self) -> bool:
        pass


class MapFlagViruses(MappingFlagBuild):
    build_name = "viruses"

    def assert_false_positive(self):

        if self.depth == 0:
            return True

        if self.depthc / self.depth > 10 and self.coverage < 5:
            return True

        return False

    def assert_vestigial(self):
        if self.mapped < 3:
            return True

        return False


class MapFlagBacteria(MappingFlagBuild):
    build_name = "bacteria"

    def assert_false_positive(self):

        if self.depth == 0:
            return True

        if self.depthc / self.depth > 10 and self.coverage < 5:
            return True

        return False

    def assert_vestigial(self):
        if self.mapped < 3:
            return True

        return False


class MapFlagProbes(MappingFlagBuild):
    build_name = "probes"

    def __init__(self, depth, depthc, coverage, mapped, windows_covered):
        super().__init__(depth, depthc, coverage, mapped, windows_covered)

        self.windows_covered = self.split_windows_covered(windows_covered)
        self.windows_covered = self.windows_covered[0] / self.windows_covered[1]

    def assert_false_positive(self):
        if self.windows_covered < 0.5:
            return True

        return False

    def assert_vestigial(self):
        if self.mapped < 3:
            return True

        return False

    @staticmethod
    def split_windows_covered(window_str: str):
        try:
            return [int(x) for x in window_str.split("/")]

        except:
            return [0, 1]
