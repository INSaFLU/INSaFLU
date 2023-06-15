
from dataclasses import dataclass

@dataclass(frozen=True)
class IntermediateFiles:
    read_classification_report: str
    contig_classification_report: str
    remap_main_report: str
    database_matches: str

    @property
    def files(self):
        return [self.read_classification_report, self.contig_classification_report, self.remap_main_report, self.database_matches]
