#!/usr/bin/env python3
import os
import json
import yaml
import anndata as ad
import numpy as np

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
METADATA_SCHEMA_DEFINITION_FILE = os.path.join(
    CURRENT_DIR, "schema_definitions", "metadata_schema_definition.yaml"
)

def get_metadata_schema_definition():
    with open(METADATA_SCHEMA_DEFINITION_FILE, "r") as f:
        schema = yaml.safe_load(f)
        sample_metadata_cols = schema["components"]["var"]["columns"]
        uns_cols = schema["components"]["uns"]["keys"]
        dataset_metadata_cols = {uns_col: "TEXT" for uns_col in uns_cols}
        dataset_metadata_attrs = {k: v["type"].upper() for k, v in schema["attributes"].items()}
        return sample_metadata_cols, dataset_metadata_cols, dataset_metadata_attrs

def get_dataset_metadata_uns(adata, uns_keys):
    dm = []
    for uns_key in uns_keys:
        v = adata.uns.get(uns_key)
        if isinstance(v, (list, np.ndarray)):
            v = ", ".join(map(str, v))
        if not isinstance(v, str):
            print(f"Warning: {uns_key} is not a string: {v}")
        dm.append(v)
    return dict(zip(uns_keys, dm))

def get_dataset_metadata_doc(adata, dataset_metadata_cols_keys, dataset_metadata_attrs, adata_file):
    dataset_metadata = get_dataset_metadata_uns(adata, dataset_metadata_cols_keys)
    dataset_metadata["protein_count"] = adata.shape[0]
    dataset_metadata["sample_count"] = adata.shape[1]
    dataset_metadata["file_name"] = os.path.basename(adata_file)
    dataset_metadata["title"] = adata.uns.get("title", "Untitled")
    return dataset_metadata

def get_sample_metadata_docs(adata, sample_metadata_cols):
    df = adata.var[sample_metadata_cols].copy()
    df["title"] = adata.uns.get("title", "Untitled")
    sample_docs = []
    for sample_id, row in df.iterrows():
        doc = row.to_dict()
        doc["sample_id"] = sample_id
        sample_docs.append(doc)
    return sample_docs

def read_ndjson(file_path):
    docs = []
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    docs.append(json.loads(line))
    return docs

def write_ndjson(file_path, docs):
    with open(file_path, "w") as f:
        for doc in docs:
            f.write(json.dumps(doc) + "\n")

def initialize_metadata_json(json_dir):
    os.makedirs(json_dir, exist_ok=True)
    dataset_file = os.path.join(json_dir, "dataset_metadata.ndjson")
    sample_file = os.path.join(json_dir, "sample_metadata.ndjson")
    open(dataset_file, "w").close()
    open(sample_file, "w").close()
    print(f"JSON files initialized in: {json_dir}")

def update_metadata_json(json_dir, adata_file, overwrite=False):
    sample_metadata_cols, dataset_metadata_cols, dataset_metadata_attrs = get_metadata_schema_definition()
    dataset_file = os.path.join(json_dir, "dataset_metadata.ndjson")
    sample_file = os.path.join(json_dir, "sample_metadata.ndjson")
    
    dataset_docs = read_ndjson(dataset_file)
    sample_docs = read_ndjson(sample_file)
    
    adata = ad.read_h5ad(adata_file)
    dataset_title = adata.uns.get("title", None)
    if not dataset_title:
        raise ValueError("No 'title' found in adata.uns")
    
    exists = any(doc.get("title") == dataset_title for doc in dataset_docs)
    if exists and not overwrite:
        raise ValueError(f"Dataset with title '{dataset_title}' already exists. Use --overwrite to update.")
    
    if exists and overwrite:
        dataset_docs = [doc for doc in dataset_docs if doc.get("title") != dataset_title]
        sample_docs = [doc for doc in sample_docs if doc.get("title") != dataset_title]
        print(f"Updating dataset with title '{dataset_title}'")
    
    new_dataset_doc = get_dataset_metadata_doc(adata, dataset_metadata_cols.keys(), dataset_metadata_attrs, adata_file)
    dataset_docs.append(new_dataset_doc)
    new_sample_docs = get_sample_metadata_docs(adata, sample_metadata_cols)
    sample_docs.extend(new_sample_docs)
    
    write_ndjson(dataset_file, dataset_docs)
    write_ndjson(sample_file, sample_docs)
    
    print(f"Metadata updated with dataset '{dataset_title}' from file '{adata_file}'.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Generate NDJSON metadata files for Sanity import from .h5ad files."
    )
    subparsers = parser.add_subparsers(dest="command")
    
    init_parser = subparsers.add_parser("initialize", help="Initialize the JSON files in a directory.")
    init_parser.add_argument("json_dir", help="Directory to store the NDJSON files.")
    
    update_parser = subparsers.add_parser("update", help="Update the JSON files with a new .h5ad file.")
    update_parser.add_argument("json_dir", help="Directory where the NDJSON files are stored.")
    update_parser.add_argument("adata_file", help="Path to the .h5ad file.")
    update_parser.add_argument("--overwrite", action="store_true", help="Overwrite if the dataset already exists.")
    
    args = parser.parse_args()
    
    if args.command == "initialize":
        initialize_metadata_json(args.json_dir)
    elif args.command == "update":
        update_metadata_json(args.json_dir, args.adata_file, args.overwrite)
    else:
        parser.print_help()
