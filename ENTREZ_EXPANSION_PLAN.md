# Expansion Plan: EntrezWrapper, Taxonomy Lineage & ORM (INSaFLU)

🛠
**Goal:** Expand `EntrezWrapper` to support **taxonomy lineage retrieval**, **organism name resolution**, and a **simple ORM + display layer** for references in INSaFLU.

Focus: **simple, install-ready schema + fast queries + cached lineage**

**Primary outcomes**

- Retrieve full lineage from NCBI TaxID
- Resolve organism names → canonical TaxIDs
- Store taxonomy in a **simple relational structure**
- Provide display-ready lineage via cached field (`lineage_path`)
</aside>

🚫
**Non-goals (for now)**

- No full graph/tree taxonomy engine
- No heavy UI work (keep minimal / functional)
- No perfect disambiguation system(iterative improvement)
- No complex biological reasoning system

---

## 1) Summary (what is being expanded)

This plan details the expansion of the `EntrezWrapper` class to:

1. **Fetch taxonomic lineages** given a TaxID
2. **Search organisms by name** using `taxonomy_utils.py` (previous class in Basespace)
3. **Design a structured ORM table** for reference storage and display
4. **Provide cached lineage** for display and filtering

---

## 2) Feedback notes (Django expansion plan)

🧩
**Key Design Decision :** **Simplified Taxonomy Model**

Instead of multiple rank tables or a custom clade system:

- We will use a **single self-referential Taxon Table** + cached lineage string in \*\*\*\*`ReferenceSource`

This avoids overengineering, while still preserving hierarchy.

### **Final chosen approach**

- One `Taxon` table (all ranks unified)
- Self-referential `parent` field
- `rank` fields stores level (kingdom, genus, species, etc..)
- `ReferenceSource` stores:
  1. FK to Taxon
  2. cached `lineage_path`

---

### ORM design (clean version)

```python
class Taxon(models.Model):
    taxid = models.IntegerField(unique=True, db_index=True)
    name = models.CharField(max_length=255, db_index=True)
    rank = models.CharField(max_length=50, db_index=True, default = TaxonConstants.NO_RANK)
    rank_raw = models.CharField(default = "no rank")

    parent = models.ForeignKey(
        "self",
        null=True,
        blank=True,
        on_delete=models.SET_NULL,
        related_name="children"
    )
```

---

```python
class ReferenceSource(models.Model):
    # existing fields...
    file = models.CharField(max_length=100, blank=True, null=True)
    description = models.CharField(max_length=300, blank=True, null=True)

    # core taxonomy link
    taxon = models.ForeignKey(Taxon, null=True, blank=True, on_delete=models.SET_NULL)

    # cached lineage (for fast display + filtering)
    lineage_path = models.TextField(blank=True)

    # optional structured cache
    lineage_json = models.JSONField(default=list, blank=True)
```

---

### Why?

- Simples schema (no rank explosion)
- Flexible for all taxonomy depths
- Compatible with NCBI hierarchy
- Fast reads via cached lineage
- No complex recursive ORM joins needed for common use cases

---

## 3) Key design decision: Lineage strategy

**Chosen approach: Hybrid storage**

We store lineage in two forms:

- **Relational truth:** `Taxon.parent`
- **Fast represetation:** `ReferenceSource.lineage_path`

- Decision recap + benefits
  - No expensive recursive queries for display
  - Still supports hierarchy traversal when needed
  - Keeps DB normalized without sacrificing performance

---

## 4) EntrezWrapper (Current Structure)

### 4.1 Existing EntrezQuery subclasses

- `EntrezFetchTaxidDescription` - get scientific name from taxid
- `EntrezFetchProteinAccession_Taxon` - map protein accessions → taxids
- `EntrezFetchAccessionDescription` - get description from accession

### 4.2 Current capabilities

- Query NCBI databases via Entrez Direct tools
- Parse tabular output into DataFrames
- Basic batching/chunking (`chunksize`)

### 4.3 What’s missing (the upgrade)

- Lineage retrieval
- Name-based searching (with fuzzy matching for malformed names)
- Biopython integration for complex queries (optional, but helpful)

---

## 5) Query classes

- 5.1 New query class: `EntrezFetchTaxidLineage`

  ```python
  class EntrezFetchTaxidLineage(EntrezQuery):
      """
      Fetch complete taxonomic lineage for a given taxid.

      Purpose: Get full hierarchy from superkingdom → species

      NCBI Database: taxonomy
      Output Fields: rank, name, taxid

      Example:
          Input taxid: 562 (E. coli)
          Output: superkingdom:Bacteria, phylum:Proteobacteria, ..., species:Escherichia coli
      """

      name = "fetch_taxid_lineage"
      db = "taxonomy"
      output_columns = ["rank", "name", "taxid"]

      def query(self, query: List[str]) -> str:
          """Build Entrez Direct command to fetch lineage"""
          # Uses efetch + xtract to parse taxonomic lineage
          # Extracts: Rank, ScientificName from LineageEx

      def read_output(self, output_path: str) -> pd.DataFrame:
          """Parse lineage output into DataFrame"""
  ```

  **Implementation notes**
  - Use `efetch -db taxonomy -format xml` (Biopython parsing is usually easier)
  - Extract `LineageEx`
  - Handle missing ranks gracefully
  - Cache results to avoid repeated NCBI calls

---

- 5.2 New query class: `EntrezSearchOrganismByName`

  ```python
  class EntrezSearchOrganismByName(EntrezQuery):
      """
      Search for organisms by scientific or common name.

      Purpose: Resolve messy organism names to canonical taxids

      NCBI Database: taxonomy
      Input: Organism name (normalized by taxonomy_utils)
      Output: taxid, scientific_name, rank

      Example:
          Input: "E. coli" → normalized → "escherichia coli"
          Output: 562, "Escherichia coli", "species"
      """

      name = "search_organism_by_name"
      db = "taxonomy"
      output_columns = ["taxid", "scientific_name", "rank"]

      def query(self, query: List[str]) -> str:
          """Build esearch command for organism name"""
          # Uses esearch with term query
          # Returns top hit (most relevant)

      def read_output(self, output_path: str) -> pd.DataFrame:
          """Parse search results"""
  ```

  **Integration with taxonomy_utils**
  - Names preprocessed by `normalize_name()` + `expand_abbreviation()`
  - Fuzzy matching fallback if exact match fails
  - Optional registry lookup

---

- 5.3 Wrapper Methods: `EntrezWrapper`

  ```python
  class EntrezWrapper:

      def fetch_lineage(self, taxids: list[str]) -> dict:
          """
          Returns structured lineage per taxid.
          """

      def search_organism_name(self, names: list[str]) -> dict:
          """
          Normalize → search → resolve taxid.
          """

      def enrich_references_dataframe(self, df):
          """
          Adds:
          - organism_name
          - taxon FK mapping
          - lineage_path (cached string)
          """
  ```

---

## 6) Data flow (simplified)

```markdown
Input file (accessions)
↓
EntrezWrapper fetches: - accession → taxid → lineage
↓
Taxon table populated (if missing)
↓
ReferenceSource created: - FK → Taxon - lineage_path cached - lineage_json optional
↓
UI reads: - taxon.name - lineage_path (fast display)
```

---

## 7) Query examples

```python
# Genus Filter
ReferenceSource.objects.filter(taxon__name="Escherichia", taxon__rank="genus")
# Species filter
ReferenceSource.objects.filter(taxon__rank="species")
# Full lineage search (fast)
ReferenceSource.objects.filter(lineage_path__icontains="Proteobacteria")
# All descendants of a taxon (if needed)
Taxon.objects.filter(parent__in=[target_taxon])
```

---

## 8) UI (basic; non-priority)

**Feedback note:** not a priority; TELEVIR References Table exists. Main concern: layout space → make lineage columns optional.

- `pathogen_identification/tables.py`
- `templates/pathogen_identification/references_table_with_lineage.html`
- `pathogen_identification/views.py`
- Features: sortable columns, search/filter, lineage display with truncation, optional columns

---

## 9) Install + migrations strategy (schema ready from install)

1. Add models to `pathogen_identification/models.py`
2. Run: `python manage.py makemigrations pathogen_identification`
3. Ensure install/setup runs: `python manage.py migrate`
4. Data population later (management command flags)

<aside>
<img src="/icons/check_green.svg" alt="/icons/check_green.svg" width="40px" />

**Result:** Tables exist immediately after setup (even if empty), so development can build on a stable schema.

</aside>

---

## 🌍 REAL EXAMPLES

### SCENARIO 1

### Registering a new reference file

Input:

```markdown
file: ecoli.fasta
accession: XBC123
```

**STEP 1 - Entrez Query**

NCBI returns:

```json
{
  "taxid": 562,
  "organism": "Escherichia coli"
}
```

**STEP 2 - Fetch lineage from Entrez**

Wrapper calls:

```python
fetch_lineage([562])
```

NCBI returns:

```json
[
  { "taxid": 2, "name": "Bacteria", "rank": "kingdom" },
  { "taxid": 1224, "name": "Proteobacteria", "rank": "phylum" },
  { "taxid": 1236, "name": "Gammaproteobacteria", "rank": "class" },
  { "taxid": 543, "name": "Enterobacteriaceae", "rank": "family" },
  { "taxid": 561, "name": "Escherichia", "rank": "genus" },
  { "taxid": 562, "name": "Escherichia coli", "rank": "species" }
]
```

**STEP 3 - Populate Taxon Table**

For each lineage item:

```python
Taxon.objects.update_or_create(
    taxid=562,
    defaults={
        "name": "Escherichia coli",
        "rank": "species",
        "parent": genus_taxon
    }
)
```

The Taxon table becomes:

```markdown
Bacteria
└── Proteobacteria
└── Gammaproteobacteria
└── Enterobacteriaceae
└── Escherichia
└── Escherichia coli
```

**STEP 4 - Create ReferenceSource**

```python
ReferenceSource.objects.create(
    file="ecoli.fasta",
    taxon=ecoli_taxon,
    lineage_path="Bacteria > Proteobacteria > ... > Escherichia coli"
)
```

### SCENARIO 2

### Another species with shared taxa

Upload:

```python
salmonella.fasta
```

NCBI returns:

```markdown
Salmonella enterica
taxid = 28901
```

**Fetch lineage:**

```markdown
Bacteria
└── Proteobacteria
└── Enterobacteriaceae
└── Salmonella
└── Salmonella enterica
```

What happens inTaxon?

```python
Taxon.objects.get_or_create(taxid=543)
```

Enterobacteriaceae is shared between E.coli and Salmonella.

✅ Reused Taxa → no lockup, only new nodes are inserted

**Resulting Tree**:

```markdown
Enterobacteriaceae
├── Escherichia
│ └── E. coli
│
└── Salmonella
└── Salmonella enterica
```

### SCENARIO 3

### search_organism_name()

User types:

```markdown
"E.coli"
```

or messy input:

```markdown
"escheria coli"
```

**STEP 1 - Normalize name**

Calls:

```python
normalize_name("E.coli")
```

Returns:

```markdown
"escherichia coli"
```

**STEP 2 - Local DB search first**

Systems tries:

```python
Taxon.objects.filter(name__iexact="Escherichia coli")
```

- **If found:**
  ```json
  {
    "taxid": 562,
    "name": "Escherichia coli"
  }
  ```
- **If not found STEP 3 - Query Entrez**
  ```python
  search_organism_name(["Escherichia coli"])
  ```
  ```json
  {
    "taxid": 562,
    "scientific_name": "Escherichia coli",
    "rank": "species"
  }
  ```
  **STEP 4 - Store Locally**
  Insert/Update Taxon:
  ```python
  Taxon.objects.update_or_create(
      taxid=562,
      defaults={...}
  )
  ```
  Return Result:
  ```json
  {
    "taxid": 562,
    "name": "Escherichia coli",
    "source": "ncbi"
  }
  ```

### SCENARIO 4

### Queries after data exists

Find all E.coli refs:

```python
ReferenceSource.objects.filter(
    taxon__name="Escherichia coli"

```

Find all refs in Enterobacteriaceae:

```python
family = Taxon.objects.get(name="Enterobacteriaceae")

ReferenceSource.objects.filter(
    taxon__in=family.get_descendants()
)
```

Display lineage quickly

```python
ref.lineage_path
```

Output:

```markdown
Bacteria > Proteobacteria > ... > Escherichia coli
```

No recursion needed.

---
