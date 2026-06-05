#!/usr/bin/env python3
"""
Update a list of UniProt IDs to their current canonical accessions.

Steps:
  1. Check each ID against the local GeneChecker for the given species.
  2. Submit invalid IDs to the UniProt ID Mapping REST API.
  3. Replace each invalid ID with its mapped counterpart (first entry if demerged).
  4. Report what changed, what was unmappable, and print the final updated list.

Usage:
    python update_uniprot_ids.py --species NCBITaxon:9606 --ids-file ids.txt
    python update_uniprot_ids.py --species homo_sapiens P12345 Q9BYF1 OLDID
"""

import sys
import os
import argparse
import time
import io

import requests
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."))

from cellxgene_schema.uniprot import GeneChecker, SupportedOrganisms

BASE = "https://rest.uniprot.org"

SPECIES_ALIASES = {
    "homo_sapiens": SupportedOrganisms.HOMO_SAPIENS,
    "human": SupportedOrganisms.HOMO_SAPIENS,
    "mus_musculus": SupportedOrganisms.MUS_MUSCULUS,
    "mouse": SupportedOrganisms.MUS_MUSCULUS,
    "rattus_norvegicus": SupportedOrganisms.RATTUS_NORVEGICUS,
    "rat": SupportedOrganisms.RATTUS_NORVEGICUS,
    "drosophila_melanogaster": SupportedOrganisms.DROSOPHILA_MELANOGASTER,
    "drosophila": SupportedOrganisms.DROSOPHILA_MELANOGASTER,
    "saccharomyces_cerevisiae": SupportedOrganisms.SACCHAROMYCES_CEREVISIAE,
    "yeast": SupportedOrganisms.SACCHAROMYCES_CEREVISIAE,
}


def resolve_species(species_str: str) -> SupportedOrganisms:
    key = species_str.lower().replace(" ", "_")
    if key in SPECIES_ALIASES:
        return SPECIES_ALIASES[key]
    for organism in SupportedOrganisms:
        if organism.value == species_str:
            return organism
    valid = list(SPECIES_ALIASES.keys()) + [o.value for o in SupportedOrganisms]
    raise ValueError(f"Unknown species '{species_str}'. Valid options: {valid}")


CHUNK_SIZE = 40


def _uniprot_map_chunk(ids: list) -> dict:
    """Run a single ID mapping job for up to CHUNK_SIZE ids. Returns {from_id: to_id}."""
    r = requests.post(f"{BASE}/idmapping/run", data={
        "from": "UniProtKB_AC-ID",
        "to": "UniProtKB",
        "ids": ",".join(ids),
    })
    r.raise_for_status()
    job_id = r.json()["jobId"]

    while True:
        s = requests.get(f"{BASE}/idmapping/status/{job_id}")
        s.raise_for_status()
        j = s.json()
        status = j.get("jobStatus")
        if status in (None, "FINISHED"):
            break
        if status == "RUNNING":
            time.sleep(1.0)
        else:
            raise RuntimeError(f"Unexpected job status: {j}")

    out = requests.get(f"{BASE}/idmapping/stream/{job_id}", params={"format": "tsv"})
    out.raise_for_status()

    df = pd.read_csv(io.StringIO(out.text), sep="\t")
    from_col, to_col = df.columns[0], df.columns[1]
    mapping = {}
    for _, row in df.iterrows():
        from_id = str(row[from_col])
        to_id = str(row[to_col])
        if from_id not in mapping:
            mapping[from_id] = to_id
    return mapping


def uniprot_map(ids: list, chunk_size: int = CHUNK_SIZE) -> dict:
    """
    Submit ids to the UniProt ID Mapping API.
    Queries longer than 40 IDs are split into chunks of 40 and results recombined.
    Returns {from_id: to_id}. For demerged entries, the first returned accession is used.
    IDs with no mapping are absent from the result.
    """
    mapping = {}
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i : i + chunk_size]
        mapping.update(_uniprot_map_chunk(chunk))
    return mapping


def update_ids(ids: list, species: SupportedOrganisms) -> tuple:
    """
    Returns (updated_ids, report) where report has keys:
      'valid'    – IDs already valid, unchanged
      'updated'  – {old: new} for IDs that were remapped
      'unmapped' – IDs the API could not resolve
    """
    checker = GeneChecker(species)

    valid = []
    invalid = []
    for uid in ids:
        (valid if checker.is_valid_id(uid) else invalid).append(uid)

    print(f"  Valid (no update needed): {len(valid)}")
    print(f"  Invalid / potentially outdated: {len(invalid)}")

    updated = {}
    unmapped = []

    if invalid:
        print(f"  Querying UniProt ID Mapping API for {len(invalid)} IDs...")
        mapping = uniprot_map(invalid)

        for uid in invalid:
            if uid in mapping:
                new_id = mapping[uid]
                if checker.is_valid_id(new_id):
                    updated[uid] = new_id
                    if new_id != uid:
                        print(f"    {uid} -> {new_id}")
                else:
                    unmapped.append(uid)
                    print(f"    {uid} -> {new_id} [NOT IN LOCAL REF]")
            else:
                unmapped.append(uid)
                print(f"    {uid} -> [NOT FOUND]")

    updated_ids = [updated.get(uid, uid) for uid in ids]
    report = {"valid": valid, "updated": updated, "unmapped": unmapped}
    return updated_ids, report


def main():
    parser = argparse.ArgumentParser(description="Update UniProt IDs to current canonical accessions")
    parser.add_argument(
        "--species", required=True,
        help="Species as NCBITaxon ID (e.g. NCBITaxon:9606) or alias (e.g. homo_sapiens, mouse)",
    )
    parser.add_argument(
        "ids", nargs="*",
        help="UniProt IDs to update (positional). Use --ids-file for large lists.",
    )
    parser.add_argument(
        "--ids-file",
        help="Plain-text file with one UniProt ID per line.",
    )
    parser.add_argument(
        "--output", "-o",
        help="Write the updated ID list to this file (one per line). Default: print to stdout.",
    )
    args = parser.parse_args()

    ids = list(args.ids)
    if args.ids_file:
        with open(args.ids_file) as f:
            ids += [line.strip() for line in f if line.strip()]

    if not ids:
        parser.error("Provide at least one UniProt ID via positional args or --ids-file.")

    species = resolve_species(args.species)
    print(f"Species: {species.name} ({species.value})")
    print(f"Total IDs: {len(ids)}\n")

    updated_ids, report = update_ids(ids, species)

    print(f"\nSummary:")
    print(f"  Unchanged : {len(report['valid'])}")
    print(f"  Updated   : {len(report['updated'])}")
    print(f"  Unmapped  : {len(report['unmapped'])}")
    if report["unmapped"]:
        print(f"  Unmapped IDs: {report['unmapped']}")

    if args.output:
        with open(args.output, "w") as f:
            f.write("\n".join(updated_ids) + "\n")
        print(f"\nUpdated IDs written to {args.output}")
    else:
        print("\nUpdated IDs:")
        for uid in updated_ids:
            print(uid)


if __name__ == "__main__":
    main()
