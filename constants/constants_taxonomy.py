class TaxonConstants:
    """
    Canonical taxonomy rank constants and abbreviation mapping.
    """

    # Canonical rank names
    RANK_DOMAIN = "domain"
    RANK_SUPERKINGDOM = "superkingdom"
    RANK_KINGDOM = "kingdom"
    RANK_SUBKINGDOM = "subkingdom"
    RANK_SUPERPHYLUM = "superphylum"
    RANK_PHYLUM = "phylum"
    RANK_SUBPHYLUM = "subphylum"
    RANK_INFRAPHYLUM = "infraphylum"
    RANK_SUPERCLASS = "superclass"
    RANK_CLASS = "class"
    RANK_SUBCLASS = "subclass"
    RANK_INFRACLASS = "infraclass"
    RANK_COHORT = "cohort"
    RANK_SUBCOHORT = "subcohort"
    RANK_SUPERORDER = "superorder"
    RANK_ORDER = "order"
    RANK_SUBORDER = "suborder"
    RANK_INFRAORDER = "infraorder"
    RANK_PARVORDER = "parvorder"
    RANK_SUPERFAMILY = "superfamily"
    RANK_FAMILY = "family"
    RANK_SUBFAMILY = "subfamily"
    RANK_TRIBE = "tribe"
    RANK_SUBTRIBE = "subtribe"
    RANK_GENUS = "genus"
    RANK_SUBGENUS = "subgenus"
    RANK_SECTION = "section"
    RANK_SUBSECTION = "subsection"
    RANK_SERIES = "series"
    RANK_SUBSERIES = "subseries"
    RANK_SPECIES_GROUP = "species group"
    RANK_SPECIES_SUBGROUP = "species subgroup"
    RANK_SPECIES = "species"
    RANK_FORMA_SPECIALIS = "forma specialis"
    RANK_SUBSPECIES = "subspecies"
    RANK_VARIETAS = "varietas"
    RANK_SUBVARIETY = "subvariety"
    RANK_FORMA = "forma"
    RANK_SEROGROUP = "serogroup"
    RANK_SEROTYPE = "serotype"
    RANK_STRAIN = "strain"
    RANK_ISOLATE = "isolate"

    NO_RANK = "no rank"


    RANK_SYNONYMS_ALLOWED = {
        "superkingdom": RANK_DOMAIN,
        "kingdom": RANK_KINGDOM,

        "phylum": RANK_PHYLUM,
        "division": RANK_PHYLUM,
        "superphylum": RANK_PHYLUM,
        "subphylum": RANK_PHYLUM,

        "class": RANK_CLASS,
        "subclass": RANK_CLASS,
        "superclass": RANK_CLASS,

        "order": RANK_ORDER,
        "superorder": RANK_ORDER,
        "suborder": RANK_ORDER,

        "family": RANK_FAMILY,
        "superfamily": RANK_FAMILY,
        "subfamily": RANK_FAMILY,

        "genus": RANK_GENUS,
        "subgenus": RANK_GENUS,

        "species": RANK_SPECIES,
    }


    # Abbreviation -> canonical name
    ABBREV_MAP = {
        "spKin": RANK_SUPERKINGDOM,
        "Kin": RANK_KINGDOM,
        "sbKin": RANK_SUBKINGDOM,
        "spPhy": RANK_SUPERPHYLUM,
        "Phy": RANK_PHYLUM,
        "sbPhy": RANK_SUBPHYLUM,
        "inPhy": RANK_INFRAPHYLUM,
        "spCla": RANK_SUPERCLASS,
        "Cla": RANK_CLASS,
        "sbCla": RANK_SUBCLASS,
        "inCla": RANK_INFRACLASS,
        "Coh": RANK_COHORT,
        "sbCoh": RANK_SUBCOHORT,
        "spOrd": RANK_SUPERORDER,
        "Ord": RANK_ORDER,
        "sbOrd": RANK_SUBORDER,
        "inOrd": RANK_INFRAORDER,
        "prOrd": RANK_PARVORDER,
        "spFam": RANK_SUPERFAMILY,
        "Fam": RANK_FAMILY,
        "sbFam": RANK_SUBFAMILY,
        "Tri": RANK_TRIBE,
        "sbTri": RANK_SUBTRIBE,
        "Gen": RANK_GENUS,
        "sbGen": RANK_SUBGENUS,
        "Sec": RANK_SECTION,
        "sbSec": RANK_SUBSECTION,
        "Ser": RANK_SERIES,
        "sbSer": RANK_SUBSERIES,
        "Sgr": RANK_SPECIES_GROUP,
        "sbSgr": RANK_SPECIES_SUBGROUP,
        "Spe": RANK_SPECIES,
        "Fsp": RANK_FORMA_SPECIALIS,
        "sbSpe": RANK_SUBSPECIES,
        "Var": RANK_VARIETAS,
        "sbVar": RANK_SUBVARIETY,
        "For": RANK_FORMA,
        "Srg": RANK_SEROGROUP,
        "Srt": RANK_SEROTYPE,
        "Str": RANK_STRAIN,
        "Iso": RANK_ISOLATE,
    }

    # Canonical name -> abbreviation
    NAME_TO_ABBREV = {v: k for k, v in ABBREV_MAP.items()}

    @classmethod
    def from_abbrev(cls, abbrev: str) -> str | None:
        """
        Convert abbreviation to canonical taxon name.
        """
        return cls.ABBREV_MAP.get(abbrev)

    @classmethod
    def to_abbrev(cls, name: str) -> str | None:
        """
        Convert canonical taxon name to abbreviation.
        """
        return cls.NAME_TO_ABBREV.get(name)

    @classmethod
    def normalize_rank(cls, value: str) -> str:
        if value is None:
            return cls.NO_RANK

        value = value.strip().lower()

        if value in cls.RANK_SYNONYMS_ALLOWED:
            return cls.RANK_SYNONYMS_ALLOWED[value]

        if value in [k.lower() for k in cls.ABBREV_MAP]:
            for abbr, rank in cls.ABBREV_MAP.items():
                if value == abbr.lower():
                    return rank

        return cls.NO_RANK

    @classmethod
    def normalize(cls, value: str) -> str | None:
        """
        Accept either a canonical name or abbreviation and return
        the canonical taxon name.
        """
        if value in cls.NAME_TO_ABBREV:
            return value

        return cls.ABBREV_MAP.get(value)

    def __iter__(self):
        """
        Iterate over canonical taxon names.
        """
        return iter(self.NAME_TO_ABBREV.keys())
    
    def __contains__(self, item: str) -> bool:
        """
        Check if a value is a valid taxon name or abbreviation.
        """
        return item in self.NAME_TO_ABBREV or item in self.ABBREV_MAP
    
    def __len__(self) -> int:
        """
        Return the number of canonical taxon ranks.
        """
        return len(self.NAME_TO_ABBREV)
    