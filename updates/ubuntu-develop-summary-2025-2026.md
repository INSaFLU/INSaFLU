# ubuntu-develop — Yearly Development Summary (May 2025 – May 2026)

**Period**: 2025-05-22 to 2026-05-22  
**Commits**: 318  
**Contributors**: SantosJGND, Joao Santos, dsobral

---

## 1. TELEVIR Pathogen Identification Module (dominant effort)

- **Pipeline refactoring**: Replaced SoftwareTree node indexing by indices with direct PK-based lookups; introduced pipeline type registration and available-node querying to avoid unnecessary build steps.
- **Aggregate Reporting**: New `FinalReportCompound` / `ReportGroup` models for grouping and aggregating classification results across runs; JSON heatmap generation/storage for clades.
- **Tree Deployment**: Significant rework — node-to-pipeline mapping, branching tree support, improved logging, safe error handling for individual leaves.
- **Reference Management**: New `ReferenceTaxid` model with lineage fields, NCBI Entrez wrapper for taxid recovery, cached reference source files, improved file filters.
- **Classifier Outputs**: Tables for tracking original classifier reports; output registration and display in intermediate reports.
- **Project Tags**: Introduced `ProjectTag` / `ProjectTagAssignment` tables and AJAX tag management.

## 2. SLURM Migration

- Adapted job tracking and wait queues for SLURM (`submit_job`, `sbatch` config).
- Removed SGE references from `process_SGE` script; moved process type definitions to constants, added job-specific memory/CPU params.
- New bash script file type for temp files.

## 3. Django & Infrastructure Upgrades

- Django 4 compatibility: replaced deprecated `is_ajax()`, updated URL patterns to `re_path`.
- Replaced `LockedAtomicTransaction` with standard `transaction.atomic` throughout.
- Migrated breadcrumbs from `bootstrap_breadcrumbs` to `view_breadcrumbs` library across all apps.
- Migrated from `python-decouple` (dotenv) to Django's native `.env` support.
- URL formatting, lazy translation updates for Django >3.2.

## 4. Metagenomics Pipeline Enhancements

- **Metaphlan integration**: Conda deployment setup, paired-end support, output parsing.
- **Voyager classifier**: PE method support, JSON taxonomy parsing, ONT technology configuration.
- **Kraken2 / Centrifuge**: Widened confidence parameter range, added unique read count to reports, phage filtering, `min-hits` defaults.
- Standardized classifier output processing across Kraken2, Centrifuge, Metaphlan, and Voyager.

## 5. Reference / Mapping / Assembly

- Added **MPXV Rivers** and **MPXV All Clades Nextstrain build** to default references.
- **RSV paired primers** table; IRMA & iVAR default settings.
- Medaka model updates, variant calling and consensus generation fixes.
- Snippy configured to use TELEVIR binary.
- BWA filter added as polyvalent software for host read filtering.
- Abricate database project updates, column count fix.
- Primer set management functionality.

## 6. UI & Frontend

- New report pages with colors and sample detail compound layout.
- Main page and project frontend edits; error rate display on report pages.
- Nextclade multiple dataset links; project settings breadcrumbs.

## 7. Refactoring & Code Quality

- Extensive linting (`lynt`) across `pathogen_identification/`, `televir/`, and utilities.
- Refactored `general_utils` into `pathutils` with explicit staticmethods.
- Removed dead code (deprecated scripts, `metagenomics_settings` pipeline step).
- Safety guards: empty taxids, missing bam files, long descriptions, empty reports.

## 8. Bug Fixes & Operations

- Fixed reads recovery in preprocessing for multi-method structure.
- Removed `delete_media` on failed pipeline node registration for safety.
- Parallel preprocessing with BWA filter.
- Database migration resets (`pathogen_identification`, `settings`).
- Kill processes on project deletion.

---

**Key architectural shift**: Migration from SGE to SLURM as the job scheduler, alongside a major refactor of the TELEVIR module to use PK-based node lookups instead of fragile index-based tree navigation, with enhanced aggregate reporting capabilities.
