## pairwise matrix by individual reads

from typing import List

from dataclasses import dataclass

from abc import ABC, abstractmethod


@dataclass
class Clade:
    """
    Clade object
    """

    name: str
    leaves: list
    private_proportion: float

    shared_proportion_std: float
    shared_proportion_min: float
    shared_proportion_max: float


class CladeFilterMethod(ABC):
    def __init__(self, reference_clade: Clade):
        self.reference_clade = reference_clade

    def filter_clades(self, clades: List[Clade]) -> List[Clade]:
        """
        Return filtered clades
        """
        return [clade for clade in clades if self.filter_clade(clade)]

    @abstractmethod
    def filter_clade(self, clade: Clade) -> bool:
        """
        Return True if clade passes filter
        """
        pass


class CladeFilterByPrivateProportion(CladeFilterMethod):
    def filter_clade(self, clade: Clade) -> bool:
        """
        Return True if clade passes filter
        """
        return clade.private_proportion >= self.reference_clade.private_proportion


class CladeFilterBySharedProportion(CladeFilterMethod):
    def filter_clade(self, clade: Clade) -> bool:
        """
        Return True if clade passes filter
        """
        return clade.shared_proportion_max >= self.reference_clade.shared_proportion_max


class CladeFilterComposed(CladeFilterMethod):
    def filter_clade(self, clade_obj: Clade) -> bool:
        if clade_obj.private_proportion < self.reference_clade.private_proportion:
            return False

        if clade_obj.shared_proportion_min < self.reference_clade.shared_proportion_min:
            return False

        return True


class CladeFilter:
    def __init__(self, reference_clade: Clade):
        self.reference_clade = reference_clade
        print("reference_clade", reference_clade)
        self.filters: List[CladeFilterMethod] = [
            CladeFilterByPrivateProportion(self.reference_clade),
            CladeFilterBySharedProportion(self.reference_clade),
        ]

    def add_filter(self, filter: CladeFilterMethod):
        """
        Add filter to filter manager
        """
        self.filters.append(filter)

    def filter_clade(self, clade: Clade) -> bool:
        """
        Return True if clade passes filter
        """
        return all([filter.filter_clade(clade) for filter in self.filters])

    def filter_clades(self, clades: List[Clade]) -> List[Clade]:
        """
        Return filtered clades
        """
        for filter in self.filters:
            clades = filter.filter_clades(clades)
        return clades
