# This script will take a h5ad file and annotate the features with the standard ontologies

import sys
import os
import argparse
import anndata as ad
import pandas as pd
import typing

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../cellxgene_schema"))
import cellxgene_schema.metadata_db as metadata_db
import cellxgene_schema.schema as schema
import cellxgene_schema.validate as validate
import env


standard_annotations = {
    "enrichment_strategy": ("Immunoprecipitation", "Immunoprecipitation"),
    "cell_type_ontology_term_id": ("CL:0000057", "fibroblast"),
    "disease_ontology_term_id": ("PATO:0000461", "normal"),
    "organism_ontology_term_id": ("NCBITaxon:9606", "Homo sapiens"),
    "sex_ontology_term_id": ("PATO:0000384", "male"),
    "development_stage_ontology_term_id": ("HsapDv:0000258", "adult"),
    "tissue_type": (pd.Categorical(["cell culture"]), "cell culture"),
}


def main():
    parser = argparse.ArgumentParser(description="Annotate a h5ad file with standard ontologies")
    parser.add_argument("h5ad_file", type=str, help="Path to the h5ad file")
    parser.add_argument("--h5ad_out", type=str, help="Path to the output h5ad file", required=False)
    args = parser.parse_args()

    adata = ad.read_h5ad(args.h5ad_file)
    validator = validate.Validator(adata)
    validator.adata = adata
    validator._set_schema_def()
    # sample_metadata_cols, dataset_metadata_cols, dataset_metadata_attrs = metadata_db.get_metadata_schema_definition()
    schema_def = schema.get_schema_definition()
    print(adata)

    ####################
    ## Annotate uns
    ####################

    if "title" in adata.uns_keys():
        print("title is in adata.uns")
    else:
        value = input("Enter the value for title: ")
        adata.uns["title"] = value

    _, dataset_metadata_cols, _ = metadata_db.get_metadata_schema_definition()
    pub_entries = [key for key in dataset_metadata_cols if key.startswith("publication")]
    pub_entries_present = [pub_entry in adata.uns for pub_entry in pub_entries]
    if all(pub_entries_present):
        print("All publication entries are present in adata.uns")
    elif any(pub_entries_present):
        for pub_entry in pub_entries[~pub_entries_present]:
            value = input(f"Enter the value for {pub_entry}: ")
            adata.uns[pub_entry] = value
    else:
        get_with_doi = input("Do you have a DOI? (y/n): ") == "y"
        if get_with_doi:
            doi = input("Enter the DOI: ")
            doi_metadata = metadata_db.get_doi_metadata(doi)
            adata.uns["publication_doi"] = doi
            for key, value in doi_metadata.items():
                adata.uns[key] = value
                print(f"Assigned {value} to {key} is in adata.uns")
        else:
            for pub_entry in pub_entries:
                value = input(f"Enter the value for {pub_entry}: ")
                adata.uns[pub_entry] = value

    ####################
    ## Annotate var
    ####################
    df = adata.var.copy()
    for col in schema_def["components"]["var"]["columns"]:
        if col in df.columns:
            print(f"{col} is in adata.var.columns")
        else:
            while True:
                value = input(
                    f"Enter the value for {col} (Default: {standard_annotations[col][0]} ({standard_annotations[col][1]})): "
                )
                if value == "":
                    value = standard_annotations[col][0]
                if isinstance(value, typing.Iterable) and len(value) == 1:
                    df[col] = value[0]
                    df[col] = df[col].astype(value.dtype)
                else:
                    df[col] = value
                #     print(col)
                column_def = validator._get_column_def("var", col)
                column = getattr(df, col)

                # First check if there are dependencies with other columns and work with a subset of the data if so
                if "dependencies" in column_def:
                    column = validator._validate_column_dependencies(df, "var", col, column_def["dependencies"])

                # If after validating dependencies there's still values in the column, validate them.
                if len(column) > 0:
                    if "warning_message" in column_def:
                        validator.warnings.append(column_def["warning_message"])
                    validator._validate_column(column, col, "var", column_def)
                if len(validator.errors) == 0:
                    adata.var[col] = df[col]
                    print(f"{col} is valid. annotated with {value}")
                    break
                else:
                    print(validator.errors)
                    print(validator.warnings)
                    validator.reset()

    write = input("Do you want to write the annotated adata to a file? (y/n): ") == "y"
    if write:
        outfile = args.h5ad_out if args.h5ad_out else input("Enter the output file path: ")
        adata.write_h5ad(outfile)


if __name__ == "__main__":
    main()
