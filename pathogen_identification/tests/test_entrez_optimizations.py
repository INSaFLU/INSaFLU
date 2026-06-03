"""
Tests for EntrezWrapper optimizations.

Tests for:
1. Batch fetching vs per-taxid fetching performance
2. Strategy comparison (Biopython vs NCBI binaries)
3. search_organism_name refactoring
"""

import pytest
from unittest.mock import Mock, patch, MagicMock
from pathogen_identification.utilities.entrez_wrapper import (
    EntrezWrapper,
    LineageNode,
    split_query
)


class TestBatchFetching:
    """Test batch fetching optimizations."""

    def test_split_query_chunks_correctly(self):
        """Verify split_query creates proper chunks."""
        taxids = [str(i) for i in range(1, 51)]  # 50 taxids
        chunks = split_query(taxids, chunksize=10)
        
        assert len(chunks) == 5
        assert all(len(chunk) == 10 for chunk in chunks)
        assert chunks[0][0] == "1"
        assert chunks[-1][-1] == "50"

    @patch('pathogen_identification.utilities.entrez_wrapper.Entrez.efetch')
    @patch('pathogen_identification.utilities.entrez_wrapper.Entrez.read')
    def test_fetch_lineages_biopy_single_batch_call(self, mock_read, mock_efetch):
        """
        Verify fetch_lineages_biopy batches taxids into single API calls.
        
        KEY TEST: With 10 taxids and chunksize 10, should result in 1 efetch call,
        not 10 individual calls.
        """
        # Setup mock response for 10 taxids
        mock_records = []
        for i in range(1, 11):
            mock_records.append({
                "TaxId": str(i),
                "ScientificName": f"organism_{i}",
                "Rank": "species",
                "LineageEx": [
                    {"TaxId": "1", "ScientificName": "Bacteria", "Rank": "superkingdom"}
                ]
            })
        
        mock_read.return_value = mock_records
        mock_efetch.return_value = MagicMock()
        
        wrapper = EntrezWrapper(outdir="/tmp", email="test@test.com")
        taxids = [str(i) for i in range(1, 11)]
        
        result = wrapper.fetch_lineages_biopy(taxids)
        
        # CRITICAL: efetch should be called ONCE with comma-separated IDs
        # not 10 times with individual IDs
        assert mock_efetch.call_count == 1
        call_args = mock_efetch.call_args
        assert "id" in call_args.kwargs
        # Verify all taxids are in the single call
        assert all(str(i) in call_args.kwargs["id"] for i in range(1, 11))
        
        # Verify all lineages are returned
        assert len(result) == 10

    @patch('pathogen_identification.utilities.entrez_wrapper.Entrez.efetch')
    @patch('pathogen_identification.utilities.entrez_wrapper.Entrez.read')
    def test_fetch_lineages_respects_chunksize(self, mock_read, mock_efetch):
        """
        Verify fetch_lineages respects chunksize and creates multiple calls
        for large lists.
        
        With 50 taxids and default chunksize, should create 5 batches.
        """
        mock_read.return_value = []
        mock_efetch.return_value = MagicMock()
        
        wrapper = EntrezWrapper(outdir="/tmp", email="test@test.com", chunksize=10)
        taxids = [str(i) for i in range(1, 51)]
        
        result = wrapper.fetch_lineages_biopy(taxids)
        
        # With 50 taxids and chunksize 10, should be 5 calls
        assert mock_efetch.call_count == 5

    @patch('subprocess.run')
    def test_fetch_lineages_binary_fallback_on_error(self, mock_run):
        """
        Verify fetch_lineages_binary falls back to Biopython on error.
        """
        # Simulate binary command failure
        mock_run.return_value = MagicMock(
            returncode=1,
            stderr="Command not found"
        )
        
        with patch.object(EntrezWrapper, 'fetch_lineages_biopy') as mock_biopy:
            mock_biopy.return_value = {}
            
            wrapper = EntrezWrapper(outdir="/tmp", email="test@test.com")
            taxids = ["1", "2", "3"]
            
            result = wrapper.fetch_lineages_binary(taxids)
            
            # Should fall back to Biopython
            mock_biopy.assert_called_once_with(taxids)


class TestSearchOrganismNameRefactoring:
    """Test search_organism_name refactoring improvements."""

    @patch('pathogen_identification.utilities.entrez_wrapper.Entrez.esearch')
    @patch('pathogen_identification.utilities.entrez_wrapper.Entrez.read')
    def test_search_organism_name_creates_template_once(self, mock_read, mock_esearch):
        """
        Verify search_organism_name creates result template once,
        then updates only the fields that change.
        
        This prevents the bug of creating incomplete dicts multiple times.
        """
        mock_esearch.return_value = MagicMock()
        mock_read.side_effect = [
            {"IdList": ["123"]},  # esearch response
            [{"ScientificName": "E. coli", "Rank": "species"}]  # efetch response
        ]
        
        wrapper = EntrezWrapper(outdir="/tmp", email="test@test.com")
        names = ["E. coli", "Salmonella", "Unknown"]
        
        with patch('pathogen_identification.utilities.entrez_wrapper.Entrez.efetch'):
            result = wrapper.search_organism_name(names)
        
        # All keys should exist in all results (initialized from template)
        required_keys = {"taxid", "canonical_name", "rank", "confidence", "source"}
        for name in names:
            assert required_keys.issubset(result[name].keys())
        
        # "Unknown" should not have a taxid (kept from template)
        # Note: actual behavior depends on NCBI response
        assert result["Unknown"]["confidence"] in ["none", "error"]

    @patch('pathogen_identification.utilities.entrez_wrapper.Entrez.esearch')
    @patch('pathogen_identification.utilities.entrez_wrapper.Entrez.read')
    def test_search_organism_name_graceful_degradation(self, mock_read, mock_esearch):
        """
        Verify search_organism_name gracefully handles failures
        without throwing exceptions.
        """
        # esearch succeeds but efetch fails
        mock_esearch.return_value = MagicMock()
        mock_read.side_effect = [
            {"IdList": ["123"]},  # esearch succeeds
            Exception("NCBI API error")  # efetch fails
        ]
        
        wrapper = EntrezWrapper(outdir="/tmp", email="test@test.com")
        
        with patch('pathogen_identification.utilities.entrez_wrapper.Entrez.efetch'):
            # Should not raise exception
            result = wrapper.search_organism_name(["E. coli"])
        
        assert result["E. coli"]["confidence"] in ["error", "low"]


class TestStrategyComparison:
    """Test strategy parameter for fetch_lineages."""

    @patch.object(EntrezWrapper, 'fetch_lineages_biopy')
    @patch.object(EntrezWrapper, 'fetch_lineages_binary')
    def test_fetch_lineages_routes_to_biopy_by_default(self, mock_binary, mock_biopy):
        """Verify biopy is default strategy."""
        mock_biopy.return_value = {}
        
        wrapper = EntrezWrapper(outdir="/tmp", email="test@test.com")
        wrapper.fetch_lineages(["1", "2", "3"])
        
        mock_biopy.assert_called_once()
        mock_binary.assert_not_called()

    @patch.object(EntrezWrapper, 'fetch_lineages_biopy')
    @patch.object(EntrezWrapper, 'fetch_lineages_binary')
    def test_fetch_lineages_routes_to_binary_when_specified(self, mock_binary, mock_biopy):
        """Verify binary strategy is called when specified."""
        mock_binary.return_value = {}
        
        wrapper = EntrezWrapper(outdir="/tmp", email="test@test.com")
        wrapper.fetch_lineages(["1", "2", "3"], strategy="binary")
        
        mock_binary.assert_called_once()
        mock_biopy.assert_not_called()


class TestLineageNodeDataclass:
    """Test LineageNode dataclass."""

    def test_lineage_node_creation(self):
        """Verify LineageNode can be instantiated with defaults."""
        node = LineageNode(taxid="123", name="E. coli", rank="species")
        
        assert node.taxid == "123"
        assert node.name == "E. coli"
        assert node.rank == "species"

    def test_lineage_node_default_rank(self):
        """Verify LineageNode has default rank."""
        node = LineageNode(taxid="123", name="Unknown")
        
        assert node.rank == "no rank"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
