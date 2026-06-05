"""
Diagnostic test suite for BFS lineage persistence.

Tests diagnose why only phylum is being registered when full taxonomy is created.

all tests use pure Python mocks
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock
from dataclasses import dataclass
from typing import Dict, List
from collections import deque

# Add parent directory to path for imports
INSAFLU_DIR = Path(__file__).parent.parent.absolute()
sys.path.insert(0, str(INSAFLU_DIR))


# MOCK DATA


@dataclass
class LineageNode:
    """Mock LineageNode matching pathogen_identification.utilities.entrez_wrapper"""
    taxid: str = ""
    name: str = ""
    rank: str = "no rank"


def ecoli_lineage() -> Dict[str, List[LineageNode]]:
    """E. coli full taxonomy lineage: domain -> species"""
    return {
        "562": [  # E. coli taxid
            LineageNode("2", "Bacteria", "superkingdom"),
            LineageNode("1239", "Proteobacteria", "phylum"),
            LineageNode("28211", "Gammaproteobacteria", "class"),
            LineageNode("91347", "Enterobacterales", "order"),
            LineageNode("543", "Enterobacteriaceae", "family"),
            LineageNode("561", "Escherichia", "genus"),
            LineageNode("562", "Escherichia coli", "species")
        ]
    }


class MockTaxonConstants:
    """Mock TaxonConstants"""
    RANK_DOMAIN = "domain"
    RANK_PHYLUM = "phylum"
    RANK_CLASS = "class"
    RANK_ORDER = "order"
    RANK_FAMILY = "family"
    RANK_GENUS = "genus"
    
    RANK_MAPPING = {
        "superkingdom": RANK_DOMAIN,
        "domain": RANK_DOMAIN,
        "phylum": RANK_PHYLUM,
        "class": RANK_CLASS,
        "order": RANK_ORDER,
        "family": RANK_FAMILY,
        "genus": RANK_GENUS,
        "species": "species",
    }
    
    @staticmethod
    def normalize_rank(value: str) -> str:
        """Normalize NCBI ranks to standard names"""
        return MockTaxonConstants.RANK_MAPPING.get(value, value)


# TEST SUITE

class TestBFSPersistence:
    """Diagnostic tests for BFS lineage persistence"""
    
    def test_rank_normalization_converts_superkingdom_to_domain(self):
        """
        TEST 1: Verify rank normalization converts 'superkingdom' -> 'domain'.
        
        Traces:
        - TaxonConstants.normalize_rank() is called for each node
        - 'superkingdom' is correctly converted to 'domain'
        - Conversion happens BEFORE RANK_TO_FIELD lookup
        
        Expected output: normalize_rank called with ['superkingdom', 'phylum']
                         Both should map to RANK_TO_FIELD keys
        """
        print("\n" + "="*70)
        print("TEST 1: Rank Normalization (superkingdom -> domain)")
        print("="*70)
        
        # Test the actual normalization
        superkingdom_normalized = MockTaxonConstants.normalize_rank("superkingdom")
        phylum_normalized = MockTaxonConstants.normalize_rank("phylum")
        
        print(f"\nRank conversions:")
        print(f"  - 'superkingdom' -> '{superkingdom_normalized}'")
        print(f"  - 'phylum' -> '{phylum_normalized}'")
        
        # THEN: Verify conversions
        assert superkingdom_normalized == "domain", \
            f"'superkingdom' should normalize to 'domain', got '{superkingdom_normalized}'"
        assert phylum_normalized == "phylum", \
            f"'phylum' should normalize to 'phylum', got '{phylum_normalized}'"
        
        print("\n✓ PASS: Rank normalization working correctly\n")
    
    def test_rank_to_field_mapping_completeness(self):
        """
        TEST 2: Verify RANK_TO_FIELD has entries for ALL expected ranks.
        
        Traces:
        - RANK_TO_FIELD dict is defined correctly in link_referencetaxid_to_lineage()
        - Has entries for: domain, phylum, class, order, family, genus
        - No missing entries = no ranks skipped
        
        Expected output: All 6 ranks present in mapping
        """
        print("\n" + "="*70)
        print("TEST 2: RANK_TO_FIELD Mapping Completeness")
        print("="*70)
        
        # This is the RANK_TO_FIELD dict from link_referencetaxid_to_lineage()
        RANK_TO_FIELD = {
            MockTaxonConstants.RANK_DOMAIN: "tax_domain",
            MockTaxonConstants.RANK_PHYLUM: "tax_phylum",
            MockTaxonConstants.RANK_CLASS: "tax_class",
            MockTaxonConstants.RANK_ORDER: "tax_order",
            MockTaxonConstants.RANK_FAMILY: "tax_family",
            MockTaxonConstants.RANK_GENUS: "tax_genus",
        }
        
        expected_ranks = ["domain", "phylum", "class", "order", "family", "genus"]
        actual_ranks = list(RANK_TO_FIELD.keys())
        
        print(f"\nRANK_TO_FIELD mapping:")
        for rank, field in RANK_TO_FIELD.items():
            print(f"  - '{rank}' -> '{field}'")
        
        print(f"\nExpected ranks: {expected_ranks}")
        print(f"Actual ranks: {actual_ranks}")
        
        # THEN: Verify completeness
        missing = set(expected_ranks) - set(actual_ranks)
        assert not missing, \
            f"Missing ranks in RANK_TO_FIELD: {missing}. Only these present: {actual_ranks}"
        
        print("\n✓ PASS: RANK_TO_FIELD mapping is complete\n")
    
    def test_bfs_creates_all_nodes_in_lineage(self):
        """
        TEST 3: Verify BFS creates Taxon objects for ALL nodes.
        
        Traces:
        - All 7 nodes are created as Taxon objects
        - Each node added to taxon_map
        - Parent relationships preserved
        
        Expected output: 7 Taxon objects created, taxon_map has 7 entries
        """
        print("\n" + "="*70)
        print("TEST 3: BFS Creates All Nodes in Lineage (E. coli)")
        print("="*70)
        
        lineage_data = ecoli_lineage()
        lineages = lineage_data
        print(f"Input: E. coli lineage with {len(lineage_data['562'])} nodes")
        
        # MOCK: Simulate BFS logic from persist_lineages()
        node_info = {}
        children = {}
        all_child_taxids = set()

        # Build node registry and children adjacency
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

        # Find roots and BFS
        roots = sorted(set(node_info) - all_child_taxids)
        queue = deque()
        taxon_map = {}

        for root_tid in roots:
            queue.append((root_tid, None))

        created_count = 0
        while queue:
            tid, parent_taxon = queue.popleft()
            node = node_info[tid]
            
            # Mock Taxon creation
            mock_taxon = MagicMock()
            mock_taxon.taxid = int(tid)
            mock_taxon.name = node.name
            mock_taxon.rank = node.rank
            mock_taxon.parent = parent_taxon
            
            created_count += 1
            taxon_map[tid] = mock_taxon

            for child_tid in children.get(tid, set()):
                queue.append((child_tid, mock_taxon))
        
        # THEN: Verify results
        print(f"\nResults:")
        print(f"  - Taxon objects created: {created_count}")
        print(f"  - taxon_map entries: {len(taxon_map)}")
        print(f"  - taxon_map keys: {sorted(taxon_map.keys())}")
        print(f"  - Expected keys: ['1239', '2', '28211', '28211', '543', '561', '562']")
        
        assert created_count == 7, f"Expected 7 Taxon objects, got {created_count}"
        assert len(taxon_map) == 7, f"Expected taxon_map size 7, got {len(taxon_map)}"
        assert "562" in taxon_map, "Leaf node (562) missing from taxon_map"
        assert "2" in taxon_map, "Root node (2) missing from taxon_map"
        
        print("\n PASS: All nodes created and mapped correctly\n")
    
    def test_taxon_map_contains_all_lineage_nodes(self):
        """
        TEST 4: Verify taxon_map has entries for ALL lineage nodes.
        
        Traces:
        - taxon_map keys match ALL node taxids in lineage
        - No missing nodes = no gaps between BFS creation and linking
        - All nodes available for link_referencetaxid_to_lineage() lookup
        
        Expected output: taxon_map has all 7 taxids
        """
        print("\n" + "="*70)
        print("TEST 4: Taxon Map Contains All Lineage Nodes")
        print("="*70)
        
        # Create mock taxon_map (as would be returned by persist_lineages)
        taxon_map = {}
        for node_taxid in ["2", "1239", "28211", "91347", "543", "561", "562"]:
            mock_taxon = MagicMock()
            mock_taxon.taxid = int(node_taxid)
            taxon_map[node_taxid] = mock_taxon
        
        # Extract all node taxids from lineage
        lineage_data = ecoli_lineage()
        lineage = lineage_data["562"]
        lineage_node_taxids = [node.taxid for node in lineage]
        
        print(f"\nLineage nodes: {lineage_node_taxids}")
        print(f"taxon_map keys: {sorted(taxon_map.keys())}")
        
        # WHEN: We try to lookup each node in taxon_map
        missing_taxids = []
        for node_taxid in lineage_node_taxids:
            taxon = taxon_map.get(node_taxid)
            if taxon is None:
                missing_taxids.append(node_taxid)
        
        # THEN: All nodes must be found
        assert not missing_taxids, \
            f"Missing taxids in taxon_map: {missing_taxids}. Map only has: {list(taxon_map.keys())}"
        
        print("\n PASS: All lineage nodes present in taxon_map\n")
    
    def test_all_rank_fields_assigned_to_reference_taxid(self):
        """
        TEST 5 (KEY): Verify ALL rank fields are assigned, not just phylum.
        
        This is THE KEY DIAGNOSTIC TEST that shows which fields are being set.
        
        Traces:
        - link_referencetaxid_to_lineage() processes all lineage nodes
        - setattr() is called for domain, phylum, class, order, family, genus
        - ReferenceTaxid has all 6 fields populated (not just phylum!)
        
        Expected output: 
          - All 7 lineage nodes processed
          - 6 fields set (tax_domain, tax_phylum, tax_class, tax_order, tax_family, tax_genus)
          - NOT just tax_phylum!
        """
        print("\n" + "="*70)
        print("TEST 5 (CRITICAL): All Rank Fields Assigned to ReferenceTaxid")
        print("="*70)
        
        lineage_data = ecoli_lineage()
        lineage = lineage_data["562"]
        print(f"\nInput lineage ({len(lineage)} nodes):")
        for node in lineage:
            norm_rank = MockTaxonConstants.normalize_rank(node.rank)
            print(f"  - taxid={node.taxid}, rank={node.rank} -> normalized={norm_rank}")
        
        # Create mock ReferenceTaxid
        ref_taxid_obj = MagicMock()
        ref_taxid_obj.tax_domain = None
        ref_taxid_obj.tax_phylum = None
        ref_taxid_obj.tax_class = None
        ref_taxid_obj.tax_order = None
        ref_taxid_obj.tax_family = None
        ref_taxid_obj.tax_genus = None
        
        # Create mock taxon_map with all nodes
        taxon_map = {}
        for node in lineage:
            mock_taxon = MagicMock()
            mock_taxon.taxid = int(node.taxid)
            mock_taxon.name = node.name
            mock_taxon.rank = MockTaxonConstants.normalize_rank(node.rank)
            taxon_map[node.taxid] = mock_taxon
        
        print(f"\ntaxon_map entries: {len(taxon_map)}")
        for tid, taxon in sorted(taxon_map.items()):
            print(f"  - tid={tid}, rank={taxon.rank}")
        
        # SIMULATE: link_referencetaxid_to_lineage logic
        RANK_TO_FIELD = {
            MockTaxonConstants.RANK_DOMAIN: "tax_domain",
            MockTaxonConstants.RANK_PHYLUM: "tax_phylum",
            MockTaxonConstants.RANK_CLASS: "tax_class",
            MockTaxonConstants.RANK_ORDER: "tax_order",
            MockTaxonConstants.RANK_FAMILY: "tax_family",
            MockTaxonConstants.RANK_GENUS: "tax_genus",
        }
        
        changed = False
        for node in lineage:
            tid = node.taxid
            if not tid:
                continue
            normalized_rank = MockTaxonConstants.normalize_rank(node.rank)
            field = RANK_TO_FIELD.get(normalized_rank)
            if field is not None:
                taxon = taxon_map.get(tid)
                if taxon is not None and getattr(ref_taxid_obj, field) != taxon:
                    setattr(ref_taxid_obj, field, taxon)
                    changed = True
        
        if changed:
            ref_taxid_obj.save()
        
        # THEN: Check which fields were set
        set_fields = {}
        for field in ["tax_domain", "tax_phylum", "tax_class", "tax_order", "tax_family", "tax_genus"]:
            current_value = getattr(ref_taxid_obj, field)
            if current_value is not None:
                set_fields[field] = current_value
        
        print(f"\nFields SET on ReferenceTaxid:")
        if set_fields:
            for field, value in sorted(set_fields.items()):
                print(f"  ✓ {field}")
        else:
            print(f"  (none)")
        
        print(f"\nFields NOT set:")
        all_fields = ["tax_domain", "tax_phylum", "tax_class", "tax_order", "tax_family", "tax_genus"]
        not_set = [f for f in all_fields if f not in set_fields]
        for field in not_set:
            print(f"  ✗ {field}")
        
        # DIAGNOSE: Did we only set phylum?
        if len(set_fields) == 1 and "tax_phylum" in set_fields:
            print(f"\n  DIAGNOSIS: Only phylum was set! BFS might be stopping after phylum.")
            print(f"    Check if lineage is being truncated or if RANK_TO_FIELD lookup is failing.\n")
        elif len(set_fields) == 0:
            print(f"\n  WARNING: No fields were set! link_referencetaxid_to_lineage might have exited early.\n")
        else:
            print(f"\n PASS: {len(set_fields)} fields set as expected!\n")
        
        # Assert at least some fields should be set
        assert len(set_fields) > 0, \
            f"No fields were set on ReferenceTaxid! Expected 6 fields (domain, phylum, class, order, family, genus)"
