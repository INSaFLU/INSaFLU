from abc import ABC, abstractmethod
from dataclasses import dataclass


class MappingFlagBuild(ABC):
    build_name = None

    @staticmethod
    @abstractmethod
    def assert_false_positive(depth, depthc, coverage, mapped):
        pass

    @staticmethod
    @abstractmethod
    def assert_vestigial(depth, depthc, coverage, mapped):
        pass


class MapFlagViruses(MappingFlagBuild):
    build_name = "viruses"

    @staticmethod
    def assert_false_positive(depth, depthc, coverage, mapped):
        if depthc / depth > 10 and coverage < 5:
            return True

        return False

    @staticmethod
    def assert_vestigial(depth, depthc, coverage, mapped):
        if mapped < 3:
            return True

        return False


class MapFlagBacteria(MappingFlagBuild):
    build_name = "bacteria"

    @staticmethod
    def assert_false_positive(depth, depthc, coverage, mapped):
        if depthc / depth > 10 and coverage < 5:
            return True

        return False

    @staticmethod
    def assert_vestigial(depth, depthc, coverage, mapped):
        if mapped < 3:
            return True

        return False


class MapFlagProbes(MappingFlagBuild):
    build_name = "probes"

    @staticmethod
    def assert_false_positive(depth, depthc, coverage, mapped):
        if depthc / depth > 10 and coverage < 5:
            return True

        return False

    @staticmethod
    def assert_vestigial(depth, depthc, coverage, mapped):
        if mapped < 3:
            return True

        return False
