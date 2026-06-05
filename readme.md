# subcellular-proteomics-curation

An adjusted version of CZI's [`cellxgene-schema`](https://github.com/chanzuckerberg/single-cell-curation),
repurposed for **subcellular proteomics** datasets. The original tool validates single-cell
*genomics* data against the cellxgene schema using GENCODE gene references; this fork swaps
gene validation for **UniProt protein** validation and trims the genomics-only machinery
(ATAC-seq, spatial/Visium, species mapping, schema migration, and the original CLI).

It is used as a **library** by the curation driver scripts in the companion
`curation_grassp` project (`validate_objects.py` and `annotate_objects.py`) — there is no
longer a `cellxgene-schema` command-line tool.

## Installation

The package lives in the `cellxgene_schema_cli/` subdirectory. Install it directly from
GitHub with pip:

```bash
pip install "git+https://github.com/czbiohub-sf/subcellular-proteomics-curation.git#subdirectory=cellxgene_schema_cli"
```

This pulls in all dependencies, including the ontology data package
[`grassp-ontology-guide`](https://github.com/czbiohub-sf/grassp-ontology-guide) (a CZ Biohub
fork of `cellxgene-ontology-guide` that adds non-animal organisms such as yeast). That wheel
is referenced directly in `requirements.txt`; if you need to install it on its own:

```bash
pip install https://github.com/czbiohub-sf/grassp-ontology-guide/releases/download/v1.9.0/cellxgene_ontology_guide-1.9.0-py3-none-any.whl
```

### Install from a local clone

```bash
make install   # runs: cd cellxgene_schema_cli && pip install -r requirements.txt && pip install .
```

## Usage

The package is consumed as a library. The key entry points:

```python
from cellxgene_schema import validate
from cellxgene_schema.write_labels import AnnDataLabelAppender
from cellxgene_schema.uniprot import GeneChecker

# Validate an in-memory AnnData object
validator = validate.Validator()
validator.adata = adata
validator._set_schema_def()
validator._deep_check()
print(validator.errors, validator.warnings)

# Append human-readable labels (protein names, lengths, locations, ontology labels)
appender = AnnDataLabelAppender(adata)
appender._add_labels()
```

See `validate_objects.py` and `annotate_objects.py` in the `curation_grassp` project for the
full validation and annotation workflow.

## Protein references

Protein validation uses gzipped UniProt TSV files bundled in
`cellxgene_schema_cli/cellxgene_schema/uniprot_files/`. Supported organisms are configured in
the `SupportedOrganisms` enum in `cellxgene_schema/uniprot.py`. The reference files are
regenerated with `cellxgene_schema_cli/scripts/protein_processing.py`, and outdated UniProt
accessions can be remapped with `scripts/update_uniprot_ids.py`.

## Schema definition

The validation rules live in
`cellxgene_schema_cli/cellxgene_schema/schema_definitions/schema_definition.yaml`.

## Development

```bash
make unit-test   # run the unit tests
make check       # pre-commit (black, ruff) + mypy
```
