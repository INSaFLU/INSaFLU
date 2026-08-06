
import sys
from unittest.mock import MagicMock

sys.modules["django"] = MagicMock()
sys.modules["django.conf"] = MagicMock()
sys.modules["django.contrib"] = MagicMock()
sys.modules["django.contrib.auth"] = MagicMock()
sys.modules["django.contrib.auth.models"] = MagicMock()
sys.modules["django.db"] = MagicMock()
sys.modules["django.db.models"] = MagicMock()

mock_models = MagicMock()
sys.modules["pathogen_identification.models"] = mock_models

from pathogen_identification.utilities.entrez_wrapper import EntrezWrapper, LineageNode


def test_full_lineage():

    wrapper = EntrezWrapper.__new__(EntrezWrapper)

    mock_ref = MagicMock()
    mock_models.ReferenceTaxid.objects.get.return_value = mock_ref

    lineage = [
        LineageNode("2", "Bacteria", "superkingdom"),
        LineageNode("1239", "Proteobacteria", "phylum"),
        LineageNode("28211", "Gammaproteobacteria", "class"),
        LineageNode("91347", "Enterobacterales", "order"),
        LineageNode("543", "Enterobacteriaceae", "family"),
        LineageNode("561", "Escherichia", "genus"),
    ]

    taxon_map = {node.taxid: MagicMock() for node in lineage}

    wrapper.link_referencetaxid_to_lineage("562", lineage, taxon_map)

    print("\nRESULT:", mock_ref.__dict__)

    assert mock_ref.tax_domain is not None
    assert mock_ref.tax_phylum is not None
    assert mock_ref.tax_class is not None
    assert mock_ref.tax_order is not None
    assert mock_ref.tax_family is not None
    assert mock_ref.tax_genus is not None
