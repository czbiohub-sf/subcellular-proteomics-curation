import sqlite3
import yaml
import urllib.request
import json
import re
from . import env
from anndata import AnnData
import anndata as ad
import os
import numpy as np


def get_metadata_schema_definition():
    with open(env.METADATA_SCHEMA_DEFINITION_FILE, "r") as f:
        schema = yaml.safe_load(f)
        sample_metadata_cols = schema["components"]["var"]["columns"]
        uns_cols = schema["components"]["uns"]["keys"]
        dataset_metadata_cols = {uns_col: "TEXT" for uns_col in uns_cols}
        dataset_metadata_attrs = {k: v["type"].upper() for k, v in schema["attributes"].items()}
        # dataset_metadata_cols.update(dataset_metadata_attrs)
        return sample_metadata_cols, dataset_metadata_cols, dataset_metadata_attrs


def _initialize_sample_metadata_table(conn: sqlite3.Connection, sample_metadata_cols: list):
    cur = conn.cursor()
    cur.execute(
        """
        CREATE TABLE sample_metadata (
            "title" TEXT,
            "sample_id" TEXT,
            {}
        )
        """.format(",".join(f'"{col}" TEXT' for col in sample_metadata_cols))
    )


def _initialize_dataset_metadata_table(conn: sqlite3.Connection, dataset_metadata_cols: dict):
    cur = conn.cursor()
    cur.execute(
        """
        CREATE TABLE dataset_metadata (
            "file_name" TEXT,
            {}
        )
        """.format(",".join(f'"{col}" TEXT' for col in dataset_metadata_cols))
    )


def initialize_metadata_db(db_file: str):
    if os.path.exists(db_file):
        # raise FileExistsError(
        #     f"""
        #     Database file '{db_file}' already exists.\n
        #     You can update the database with the following command:\n\n
        #     cellxgene_schema update-metadata-db {db_file} <h5ad_file> \n\n
        #     If you want to make a new database, please delete the existing one or choose a different name.
        #     """
        # )
        os.remove(db_file)
    sample_metadata_cols, dataset_metadata_cols, dataset_metadata_attrs = get_metadata_schema_definition()
    dataset_metadata_cols.update(dataset_metadata_attrs)
    with sqlite3.connect(db_file) as con:
        _initialize_sample_metadata_table(con, sample_metadata_cols)
        _initialize_dataset_metadata_table(con, dataset_metadata_cols)


def _initialize_dataset_metadata_table(conn: sqlite3.Connection, dataset_metadata_cols: dict):
    cur = conn.cursor()
    cur.execute(
        """
        CREATE TABLE dataset_metadata (
            {}
            )
            """.format(",".join(f'"{col}" TEXT' for col in dataset_metadata_cols))
    )


def check_if_sample_metadata_table_exists(conn: sqlite3.Connection):
    cur = conn.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='sample_metadata'")
    return cur.fetchone() is not None


def check_if_dataset_metadata_table_exists(conn: sqlite3.Connection):
    cur = conn.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='dataset_metadata'")
    return cur.fetchone() is not None


def update_metadata_db(db_file: str, adata_file: str, overwrite: bool = False):
    sample_metadata_cols, dataset_metadata_cols, dataset_metadata_attrs = get_metadata_schema_definition()

    if not os.path.exists(db_file):
        raise FileNotFoundError(
            f"Database file '{db_file}' does not exist. You can create a new database with the following command:\n\ncellxgene_schema initialize-metadata-db {db_file}"
        )
    with sqlite3.connect(db_file, timeout=30, isolation_level=None) as conn:
        adata = ad.read_h5ad(adata_file)
        if check_if_title_exists(conn, adata.uns["title"]):
            if not overwrite:
                raise ValueError(
                    f"Title '{adata.uns['title']}' already exists in the database. If you want to overwrite it, please use the --overwrite flag."
                )
            else:
                print(f"Warning: Overwriting title '{adata.uns['title']}' in the database.")
                cur = conn.cursor()
                cur.execute("DELETE FROM dataset_metadata WHERE title = ?", (adata.uns["title"],))
                cur.execute("DELETE FROM sample_metadata WHERE title = ?", (adata.uns["title"],))
        add_sample_metadata(conn, adata, sample_metadata_cols)
        add_dataset_metadata(conn, adata, dataset_metadata_cols, dataset_metadata_attrs, adata_file)


def add_sample_metadata(conn: sqlite3.Connection, adata: AnnData, sample_metadata_cols: list):
    sample_metadata = adata.var[sample_metadata_cols].copy()
    sample_metadata.loc[:, "title"] = adata.uns["title"]

    # Create table if it doesn't exist
    if not check_if_sample_metadata_table_exists(conn):
        _initialize_sample_metadata_table(conn, sample_metadata.columns)

    sample_metadata.to_sql(
        "sample_metadata", conn, if_exists="append", index=True, index_label="sample_id", method="multi"
    )


def add_dataset_metadata(
    conn: sqlite3.Connection, adata: AnnData, dataset_metadata_cols: dict, dataset_metadata_attrs: dict, adata_file: str
):
    dataset_metadata = get_dataset_metadata_uns(adata, dataset_metadata_cols.keys())
    dataset_metadata["protein_count"] = adata.shape[0]
    dataset_metadata["sample_count"] = adata.shape[1]
    dataset_metadata["file_name"] = os.path.basename(adata_file)
    assert set(dataset_metadata.keys()) == set(dataset_metadata_cols.keys()).union(set(dataset_metadata_attrs.keys()))

    if not check_if_dataset_metadata_table_exists(conn):
        _initialize_dataset_metadata_table(conn, dataset_metadata_cols)

    cur = conn.cursor()
    columns = ", ".join(dataset_metadata.keys())
    placeholders = ", ".join("?" * len(dataset_metadata))
    sql = f"INSERT INTO dataset_metadata ({columns}) VALUES ({placeholders})"
    cur.execute(sql, list(dataset_metadata.values()))


def get_dataset_metadata_uns(adata, uns_keys):
    dm = []
    for uns_key in uns_keys:
        v = adata.uns[uns_key]
        if isinstance(v, (list, np.ndarray)):
            v = ", ".join(v)
        if not isinstance(v, str):
            print(f"Warning: {uns_key} is not a string: {v}")
        dm.append(v)
    return dict(zip(uns_keys, dm))


def check_if_title_exists(conn: sqlite3.Connection, title: str):
    cur = conn.cursor()
    cur.execute("SELECT COUNT(*) FROM dataset_metadata WHERE title = ?", (title,))
    return cur.fetchone()[0] > 0


def delete_metadata_db_entry(db_file: str, title: str):
    with sqlite3.connect(db_file, timeout=30, isolation_level=None) as con:
        cur = con.cursor()

        # Check if the title exists
        cur.execute("SELECT COUNT(*) FROM dataset_metadata WHERE title = ?", (title,))
        if cur.fetchone()[0] == 0:
            raise ValueError(f"Title '{title}' does not exist in the database.")

        # Delete the entry
        cur.execute("DELETE FROM dataset_metadata WHERE title = ?", (title,))
        cur.execute("DELETE FROM sample_metadata WHERE title = ?", (title,))


def get_doi_metadata(doi_input):
    """
    Fetch metadata (date, authors, title, journal name) for a given DOI or DOI URL using CrossRef API.

    Args:
        doi_input (str): The DOI or DOI URL of the article.

    Returns:
        dict: A dictionary containing the publication date, authors, title, and journal name.
    """
    # Extract DOI from a URL if necessary
    doi_pattern = r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]+"
    match = re.search(doi_pattern, doi_input)
    if not match:
        return {"error": "Invalid DOI format"}

    doi = match.group(0)
    url = f"https://api.crossref.org/works/{doi}"

    try:
        with urllib.request.urlopen(url) as response:
            data = json.loads(response.read().decode())

        message = data.get("message", {})

        # Extract metadata
        title = message.get("title", ["Unknown Title"])[0]
        journal = message.get("container-title", ["Unknown Journal"])[0]
        date_parts = message.get("issued", {}).get("date-parts", [[None]])
        publication_date = "-".join(map(str, date_parts[0])) if date_parts[0] else "Unknown Date"

        # Extract authors
        authors = message.get("author", [])
        author_list = [f"{author.get('given', '')} {author.get('family', '')}".strip() for author in authors]

        return {
            "publication_title": title,
            "publication_journal": journal,
            "publication_date": publication_date,
            "publication_authors": author_list,
        }

    except Exception as e:
        return {"error": str(e)}
