from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class IntermediateFiles:
    read_classification_report: Optional[str]
    contig_classification_report: Optional[str]
    remap_main_report: Optional[str]
    database_matches: str

    @property
    def files(self):

        files = [
            self.read_classification_report,
            self.contig_classification_report,
            self.remap_main_report,
            self.database_matches,
        ]
        return [f for f in files if f is not None]
