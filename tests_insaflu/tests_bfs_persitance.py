"""
Diagnostic Test Suite: BFS Lineage Persistence

PROBLEM: Taxon objects are created correctly, but ReferenceTaxid is
only registering tax_phylum. All other rank fields (tax_domain,
tax_class, tax_order, tax_family, tax_genus) are NOT being set.

Mock functions EXACTLY mirror the real code paths:
  - simulate_persist_lineages        → persist_lineages()
  - simulate_link_referencetaxid_to_lineage → link_referencetaxid_to_lineage()

Differences from real code:
  - Uses MockTaxon/MockReferenceTaxid instead of Django models
  - No DB reads/writes
  - Returns fields_set list instead of saving to DB
"""

import pytest
from unittest.mock import Mock, MagicMock
from dataclasses import dataclass
from collections import defaultdict, deque


# Mock TaxonConstants (mirrors constants/constants_taxonomy.py) 

class TaxonConstants:
    RANK_DOMAIN = "domain"
    RANK_KINGDOM = "kingdom"
    RANK_PHYLUM = "phylum"
    RANK_CLASS = "class"
    RANK_ORDER = "order"
    RANK_FAMILY = "family"
    RANK_GENUS = "genus"
    RANK_SPECIES = "species"
    NO_RANK = "no rank"

    RANK_SYNONYMS = {
        "domain": RANK_DOMAIN,
        "superkingdom": RANK_DOMAIN,
        "kingdom": RANK_KINGDOM,
        "phylum": RANK_PHYLUM,
        "division": RANK_PHYLUM,
        "class": RANK_CLASS,
        "order": RANK_ORDER,
        "family": RANK_FAMILY,
        "genus": RANK_GENUS,
        "species": RANK_SPECIES,
        "no rank": NO_RANK,
        "clade": NO_RANK,
    }

    @staticmethod
    def normalize_rank(rank_str):
        if rank_str is None:
            return TaxonConstants.NO_RANK
        normalized = TaxonConstants.RANK_SYNONYMS.get(
            rank_str.lower(), TaxonConstants.NO_RANK
        )
        return normalized


# Field mapping (mirrors the RANK_TO_FIELD inside link_referencetaxid_to_lineage)

RANK_TO_FIELD = {
    TaxonConstants.RANK_DOMAIN: "tax_domain",
    TaxonConstants.RANK_PHYLUM: "tax_phylum",
    TaxonConstants.RANK_CLASS: "tax_class",
    TaxonConstants.RANK_ORDER: "tax_order",
    TaxonConstants.RANK_FAMILY: "tax_family",
    TaxonConstants.RANK_GENUS: "tax_genus",
}

ALL_RANK_FIELDS = [
    "tax_domain", "tax_phylum", "tax_class",
    "tax_order", "tax_family", "tax_genus",
]


@dataclass
class LineageNode:
    """Mirrors dataclass in entrez_wrapper.py"""
    taxid: str = ""
    name: str = ""
    rank: str = "no rank"


class MockTaxon:
    """Mirrors Django Taxon model"""
    def __init__(self, taxid, name, rank, parent=None):
        self.taxid = taxid
        self.name = name
        self.rank = rank
        self.parent = parent

    def __repr__(self):
        return f"Taxon(taxid={self.taxid}, rank='{self.rank}')"


class MockReferenceTaxid:
    """Mirrors Django ReferenceTaxid model"""
    def __init__(self, **kwargs):
        for field in ALL_RANK_FIELDS:
            setattr(self, field, kwargs.get(field, None))

    def __repr__(self):
        sets = {f: getattr(self, f) for f in ALL_RANK_FIELDS}
        return f"ReferenceTaxid({sets})"



def simulate_persist_lineages(lineages):
    """
    EXACT mirror of entrez_wrapper.py persist_lineages().

    Args:
        lineages: Dict[str, List[LineageNode]] — leaf taxid -> lineage list

    Returns:
        taxon_map: Dict[str, MockTaxon] — tid string -> MockTaxon
    """
    node_info = {}
    children = {}
    all_child_taxids = set()

    for leaf_taxid, lineage_list in lineages.items():
        prev = None
        for node in lineage_list:
            tid = node.taxid
            if not tid:
                continue
            if tid not in node_info:
                node_info[tid] = node
            if prev is not None:
                children.setdefault(prev, set()).add(tid)
                all_child_taxids.add(tid)
            prev = tid

    roots = sorted(set(node_info) - all_child_taxids)
    queue = deque()
    taxon_map = {}

    for root_tid in roots:
        queue.append((root_tid, None))

    while queue:
        tid, parent_taxon = queue.popleft()
        node = node_info[tid]

        taxon_obj = MockTaxon(
            taxid=int(tid),
            name=node.name,
            rank=node.rank,
            parent=parent_taxon,
        )
        taxon_map[tid] = taxon_obj

        for child_tid in children.get(tid, set()):
            queue.append((child_tid, taxon_obj))

    return taxon_map


def simulate_link_referencetaxid_to_lineage(
    ref_taxid_obj, lineage, taxon_map):
    """
    EXACT mirror of entrez_wrapper.py link_referencetaxid_to_lineage().

    Args:
        ref_taxid_obj: MockReferenceTaxid instance
        lineage: List[LineageNode] — specific lineage for this taxid
        taxon_map: Dict[str, MockTaxon] — from simulate_persist_lineages()

    Returns:
        fields_set: List of field names that were assigned
    """
    fields_set = []
    for node in lineage:
        tid = node.taxid
        if not tid:
            continue
        normalized_rank = TaxonConstants.normalize_rank(node.rank)
        field = RANK_TO_FIELD.get(normalized_rank)
        if field is not None:
            taxon = taxon_map.get(tid)
            if taxon is not None and getattr(ref_taxid_obj, field) != taxon:
                setattr(ref_taxid_obj, field, taxon)
                fields_set.append(field)
    return fields_set


def run_full_pipeline(lineages):
    """
    Run persist + link for all lineages in a single call.

    Returns dict of leaf_taxid -> fields_set for each lineage.
    """
    taxon_map = simulate_persist_lineages(lineages)
    results = {}
    for leaf_tid, lineage_list in lineages.items():
        ref = MockReferenceTaxid()
        fields = simulate_link_referencetaxid_to_lineage(
            ref, lineage_list, taxon_map
        )
        results[leaf_tid] = {
            "ref": ref,
            "fields_set": fields,
            "taxon_map": taxon_map,
        }
    return results


def ecoli_lineage():
    """Standard E. coli full lineage (all 6 ranks present)."""
    return [
        LineageNode(taxid="2", name="Bacteria", rank="superkingdom"),
        LineageNode(taxid="1239", name="Proteobacteria", rank="phylum"),
        LineageNode(taxid="28211", name="Gammaproteobacteria", rank="class"),
        LineageNode(taxid="91347", name="Enterobacterales", rank="order"),
        LineageNode(taxid="543", name="Enterobacteriaceae", rank="family"),
        LineageNode(taxid="561", name="Escherichia", rank="genus"),
        LineageNode(taxid="562", name="Escherichia coli", rank="species"),
    ]


def virus_lineage():
    """Virus lineage with only domain/family/genus (typical)."""
    return [
        LineageNode(taxid="10239", name="Viruses", rank="superkingdom"),
        LineageNode(taxid="11118", name="Flaviviridae", rank="family"),
        LineageNode(taxid="11051", name="Flavivirus", rank="genus"),
        LineageNode(taxid="12637", name="Dengue virus", rank="species"),
    ]


def lineage_with_cellular_organisms():
    """NCBI-style lineage with 'cellular organisms' no-rank root."""
    return [
        LineageNode(taxid="131567", name="cellular organisms", rank="no rank"),
        LineageNode(taxid="2", name="Bacteria", rank="superkingdom"),
        LineageNode(taxid="1239", name="Proteobacteria", rank="phylum"),
        LineageNode(taxid="28211", name="Gammaproteobacteria", rank="class"),
        LineageNode(taxid="91347", name="Enterobacterales", rank="order"),
        LineageNode(taxid="543", name="Enterobacteriaceae", rank="family"),
        LineageNode(taxid="561", name="Escherichia", rank="genus"),
        LineageNode(taxid="562", name="Escherichia coli", rank="species"),
    ]


def lineage_with_empty_taxid_in_middle():
    """Domain node OK, but an intermediate node has empty taxid."""
    return [
        LineageNode(taxid="2", name="Bacteria", rank="superkingdom"),
        LineageNode(taxid="", name="", rank=""),  # ← empty taxid
        LineageNode(taxid="1239", name="Proteobacteria", rank="phylum"),
        LineageNode(taxid="28211", name="Gammaproteobacteria", rank="class"),
        LineageNode(taxid="91347", name="Enterobacterales", rank="order"),
        LineageNode(taxid="543", name="Enterobacteriaceae", rank="family"),
        LineageNode(taxid="561", name="Escherichia", rank="genus"),
        LineageNode(taxid="562", name="Escherichia coli", rank="species"),
    ]


def lineage_with_empty_taxid_at_root():
    """The root/domain node has an empty taxid."""
    return [
        LineageNode(taxid="", name="", rank=""),  # ← empty taxid
        LineageNode(taxid="1239", name="Proteobacteria", rank="phylum"),
        LineageNode(taxid="28211", name="Gammaproteobacteria", rank="class"),
        LineageNode(taxid="91347", name="Enterobacterales", rank="order"),
        LineageNode(taxid="543", name="Enterobacteriaceae", rank="family"),
        LineageNode(taxid="561", name="Escherichia", rank="genus"),
        LineageNode(taxid="562", name="Escherichia coli", rank="species"),
    ]


#  TESTS

class TestRankNormalization:
    """Rank normalization and mapping are correct (foundational)."""

    def test_superkingdom_to_domain(self):
        result = TaxonConstants.normalize_rank("superkingdom")
        assert result == TaxonConstants.RANK_DOMAIN

    def test_phylum_stays_phylum(self):
        assert TaxonConstants.normalize_rank("phylum") == "phylum"

    def test_class_stays_class(self):
        assert TaxonConstants.normalize_rank("class") == "class"

    def test_order_stays_order(self):
        assert TaxonConstants.normalize_rank("order") == "order"

    def test_family_stays_family(self):
        assert TaxonConstants.normalize_rank("family") == "family"

    def test_genus_stays_genus(self):
        assert TaxonConstants.normalize_rank("genus") == "genus"

    def test_no_rank_stays_no_rank(self):
        assert TaxonConstants.normalize_rank("no rank") == "no rank"

    def test_none_returns_no_rank(self):
        assert TaxonConstants.normalize_rank(None) == "no rank"

    def test_unknown_rank_returns_no_rank(self):
        assert TaxonConstants.normalize_rank("bogus") == "no rank"


class TestRANK_TO_FIELD:
    """RANK_TO_FIELD has all 6 needed mappings."""

    def test_all_mappings_present(self):
        required = {
            TaxonConstants.RANK_DOMAIN: "tax_domain",
            TaxonConstants.RANK_PHYLUM: "tax_phylum",
            TaxonConstants.RANK_CLASS: "tax_class",
            TaxonConstants.RANK_ORDER: "tax_order",
            TaxonConstants.RANK_FAMILY: "tax_family",
            TaxonConstants.RANK_GENUS: "tax_genus",
        }
        for rank, expected_field in required.items():
            assert rank in RANK_TO_FIELD, f"Missing rank: {rank}"
            assert RANK_TO_FIELD[rank] == expected_field, \
                f"{rank} → {RANK_TO_FIELD[rank]}, expected {expected_field}"


class TestPersistLineages:
    """simulate_persist_lineages (mirrors real persist_lineages())."""

    def test_creates_all_nodes(self):
        lineages = {"562": ecoli_lineage()}
        taxon_map = simulate_persist_lineages(lineages)
        assert len(taxon_map) == 7  # 7 nodes for E. coli

    def test_rank_stored_raw_not_normalized(self):
        lineages = {"562": ecoli_lineage()}
        taxon_map = simulate_persist_lineages(lineages)
        # Real code stores node.rank as-is (no normalization)
        assert taxon_map["2"].rank == "superkingdom"

    def test_empty_taxid_skipped_from_map(self):
        lineages = {"562": lineage_with_empty_taxid_in_middle()}
        taxon_map = simulate_persist_lineages(lineages)
        assert "" not in taxon_map  # empty taxid never stored

    def test_all_nodes_except_empty_present(self):
        lineages = {"562": lineage_with_empty_taxid_in_middle()}
        taxon_map = simulate_persist_lineages(lineages)
        expected = {"2", "1239", "28211", "91347", "543", "561", "562"}
        for tid in expected:
            assert tid in taxon_map, f"Missing taxid {tid} in taxon_map"

    def test_empty_root_reattaches_children(self):
        lineages = {"562": lineage_with_empty_taxid_at_root()}
        taxon_map = simulate_persist_lineages(lineages)
        # root node empty -> not in taxon_map
        assert "" not in taxon_map
        # but all subsequent nodes should be present
        expected = {"1239", "28211", "91347", "543", "561", "562"}
        for tid in expected:
            assert tid in taxon_map, f"Missing taxid {tid} in taxon_map"

    def test_no_rank_node_present_in_map(self):
        lineages = {"562": lineage_with_cellular_organisms()}
        taxon_map = simulate_persist_lineages(lineages)
        assert "131567" in taxon_map

    def test_shared_ancestors_deduplicated(self):
        # Two E. coli strains sharing the same lineage
        lineages = {
            "562": ecoli_lineage(),
            "316385": [
                LineageNode(taxid="2", name="Bacteria", rank="superkingdom"),
                LineageNode(taxid="1239", name="Proteobacteria", rank="phylum"),
                LineageNode(taxid="28211", name="Gammaproteobacteria", rank="class"),
                LineageNode(taxid="91347", name="Enterobacterales", rank="order"),
                LineageNode(taxid="543", name="Enterobacteriaceae", rank="family"),
                LineageNode(taxid="561", name="Escherichia", rank="genus"),
                LineageNode(taxid="316385", name="Escherichia coli O157:H7", rank="species"),
            ],
        }
        taxon_map = simulate_persist_lineages(lineages)
        assert len(taxon_map) == 8  # 7 shared + 1 extra leaf
        assert "562" in taxon_map
        assert "316385" in taxon_map

    def test_bfs_adjacency_maintains_parent_chain(self):
        lineages = {"562": ecoli_lineage()}
        taxon_map = simulate_persist_lineages(lineages)
        # BFS starts from root, so leaf should have proper ancestors
        leaf = taxon_map["562"]
        assert leaf.parent is not None
        assert leaf.parent.taxid == 561


class TestLinkReferenceTaxid:
    """simulate_link_referencetaxid_to_lineage (mirrors real method)."""

    def test_all_six_fields_set_for_complete_lineage(self):
        lineages = {"562": ecoli_lineage()}
        taxon_map = simulate_persist_lineages(lineages)
        ref = MockReferenceTaxid()
        fields_set = simulate_link_referencetaxid_to_lineage(
            ref, lineages["562"], taxon_map
        )
        assert len(fields_set) == 6, \
            f"Expected 6 fields, got {len(fields_set)}: {fields_set}"
        for field in ALL_RANK_FIELDS:
            assert getattr(ref, field) is not None, f"{field} was not set"

    def test_virus_sets_only_available_ranks(self):
        lineages = {"12637": virus_lineage()}
        taxon_map = simulate_persist_lineages(lineages)
        ref = MockReferenceTaxid()
        fields_set = simulate_link_referencetaxid_to_lineage(
            ref, lineages["12637"], taxon_map
        )
        # Virus: superkingdom → family → genus → species
        # Should set: tax_domain, tax_family, tax_genus
        assert "tax_domain" in fields_set
        assert "tax_family" in fields_set
        assert "tax_genus" in fields_set
        # Should NOT set: phylum, class, order
        assert "tax_phylum" not in fields_set
        assert "tax_class" not in fields_set
        assert "tax_order" not in fields_set
        assert len(fields_set) == 3

    def test_cellular_organisms_root_still_sets_all_fields(self):
        lineages = {"562": lineage_with_cellular_organisms()}
        taxon_map = simulate_persist_lineages(lineages)
        ref = MockReferenceTaxid()
        fields_set = simulate_link_referencetaxid_to_lineage(
            ref, lineages["562"], taxon_map
        )
        assert len(fields_set) == 6, \
            f"Expected 6 fields, got {len(fields_set)}: {fields_set}"

    def test_empty_taxid_in_middle_reconnects_chain(self):
        lineages = {"562": lineage_with_empty_taxid_in_middle()}
        taxon_map = simulate_persist_lineages(lineages)
        ref = MockReferenceTaxid()
        fields_set = simulate_link_referencetaxid_to_lineage(
            ref, lineages["562"], taxon_map
        )
        # The empty taxid node is between "Bacteria" (tax_domain) and
        # "Proteobacteria" (tax_phylum). BFS reconnects adjacency, so
        # all 6 ranked nodes should be reachable → all 6 fields set.
        assert len(fields_set) == 6, \
            f"Expected 6 fields, got {len(fields_set)}: {fields_set}"

    def test_empty_root_skips_domain_field(self):
        lineages = {"562": lineage_with_empty_taxid_at_root()}
        taxon_map = simulate_persist_lineages(lineages)
        ref = MockReferenceTaxid()
        fields_set = simulate_link_referencetaxid_to_lineage(
            ref, lineages["562"], taxon_map
        )
        # Domain/superkingdom node had empty taxid → not in taxon_map
        # → tax_domain cannot be set. But phylum/class/order/family/genus
        # are all present and set.
        assert "tax_domain" not in fields_set, \
            "tax_domain should NOT be set (empty root taxid was never persisted)"
        assert "tax_phylum" in fields_set
        assert "tax_class" in fields_set
        assert "tax_order" in fields_set
        assert "tax_family" in fields_set
        assert "tax_genus" in fields_set
        assert len(fields_set) == 5, \
            f"Expected 5 fields, got {len(fields_set)}: {fields_set}"

    def test_node_missing_from_taxon_map_silently_skipped(self):
        # Simulate: taxon_map deliberately missing a mid-lineage taxid
        lineages = {"562": ecoli_lineage()}
        taxon_map = simulate_persist_lineages(lineages)
        # Remove class node from taxon_map to simulate BFS unreachability
        del taxon_map["28211"]
        ref = MockReferenceTaxid()
        fields_set = simulate_link_referencetaxid_to_lineage(
            ref, lineages["562"], taxon_map
        )
        assert "tax_class" not in fields_set
        assert "tax_domain" in fields_set
        assert "tax_phylum" in fields_set
        assert "tax_order" in fields_set
        assert "tax_family" in fields_set
        assert "tax_genus" in fields_set
        assert len(fields_set) == 5

    def test_field_not_overwritten_when_same_object(self):
        lineages = {"562": ecoli_lineage()}
        taxon_map = simulate_persist_lineages(lineages)
        # Pre-set one field to the same object it would get
        ref = MockReferenceTaxid(
            tax_domain=taxon_map["2"]  # same object that would be set
        )
        fields_set = simulate_link_referencetaxid_to_lineage(
            ref, lineages["562"], taxon_map
        )
        # getattr guard: if existing field == taxon → skip
        assert "tax_domain" not in fields_set, \
            "tax_domain should NOT be in fields_set (already same object)"
        assert len(fields_set) == 5


class TestFullPipeline:
    """End-to-end: persist then link for realistic scenarios."""

    def test_ecoli_e2e_all_six(self):
        results = run_full_pipeline({"562": ecoli_lineage()})
        fields = results["562"]["fields_set"]
        assert len(fields) == 6

    def test_virus_e2e_three(self):
        results = run_full_pipeline({"12637": virus_lineage()})
        fields = results["12637"]["fields_set"]
        assert len(fields) == 3

    def test_cellular_organisms_e2e_all_six(self):
        results = run_full_pipeline({"562": lineage_with_cellular_organisms()})
        fields = results["562"]["fields_set"]
        assert len(fields) == 6

    def test_empty_middle_e2e_all_six(self):
        results = run_full_pipeline({"562": lineage_with_empty_taxid_in_middle()})
        fields = results["562"]["fields_set"]
        assert len(fields) == 6

    def test_empty_root_e2e_five_only(self):
        results = run_full_pipeline({"562": lineage_with_empty_taxid_at_root()})
        fields = results["562"]["fields_set"]
        assert len(fields) == 5
        assert "tax_domain" not in fields

    def test_two_organisms_shared_ancestors(self):
        strain2_lineage = [
            LineageNode(taxid="2", name="Bacteria", rank="superkingdom"),
            LineageNode(taxid="1239", name="Proteobacteria", rank="phylum"),
            LineageNode(taxid="28211", name="Gammaproteobacteria", rank="class"),
            LineageNode(taxid="91347", name="Enterobacterales", rank="order"),
            LineageNode(taxid="543", name="Enterobacteriaceae", rank="family"),
            LineageNode(taxid="561", name="Escherichia", rank="genus"),
            LineageNode(taxid="316385", name="Escherichia coli O157:H7", rank="species"),
        ]
        results = run_full_pipeline({
            "562": ecoli_lineage(),
            "316385": strain2_lineage,
        })
        assert len(results["562"]["fields_set"]) == 6
        assert len(results["316385"]["fields_set"]) == 6
        # Both should share the same ancestor taxon objects
        assert results["562"]["ref"].tax_domain is results["316385"]["ref"].tax_domain
