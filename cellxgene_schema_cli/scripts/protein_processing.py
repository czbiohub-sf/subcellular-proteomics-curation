import requests
import csv
import gzip

# import hashlib
import os

# import statistics
import sys

# import urllib.request
from collections import defaultdict
from typing import Dict

# import gtf_tools
import yaml

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../cellxgene_schema"))
import env

import re
import requests
from requests.adapters import HTTPAdapter, Retry


class ProteinProcessingResult:
    def __init__(
        self, protein_id: str, protein_name: str, protein_version: str, protein_length: str, protein_type: str
    ):
        self.protein_id = protein_id  # ex: ENSG00000141510, should be unique
        self.protein_name = protein_name  # ex: TP53, not necessarily unique
        self.protein_version = protein_version
        self.protein_length = protein_length
        self.protein_type = protein_type


class ProteinProcessor:
    def __init__(self):
        # Global mapping of protein id to all the metadata we need to proteinrate output files
        self.protein_metadata: dict[str, ProteinProcessingResult] = {}
        # Mapping from description (ex: "mus musculus") to a list of protein ids
        self.protein_ids_by_description: dict[str, list[str]] = {}

    # def write_gzip(self, data: str, output_filename: str):
    #     """
    #     Writes data to a gziped file. The date modified is not written to the gzip file. This allows
    #     comparing the hashed contents of two gziped files to determine if they are the same.

    #     :param str data: data to write
    #     :param str output_file: path to output file

    #     :rtype: None
    #     """
    #     with open(output_filename, "wb") as fileobj, gzip.GzipFile(mode="wb", fileobj=fileobj, mtime=0) as myzip:
    #         myzip.write(data.encode("utf-8"))

    # def digest(self, file_name: str) -> str:
    #     with open(file_name, "rb") as f:
    #         # Read the contents of the file in chunks
    #         chunk_size = 1024
    #         hasher = hashlib.sha256()
    #         while chunk := f.read(chunk_size):
    #             hasher.update(chunk)
    #     return hasher.hexdigest()

    # def _parse_gtf(self, gtf_path: str, protein_info_description: str):
    #     """
    #     Parses a gziped GTF file to get protein and transcript info into a gziped comma-separated file with the following
    #     structure, with three columns and no header: 1) protein id, 2) protein,
    #     3) protein

    #     :param str gtf_path: path to gzipped gtf file
    #     :param str output_json_file: path to output json

    #     :rtype: None
    #     """

    #     protein_lengths = self._get_protein_lengths_from_gtf(gtf_path)
    #     with gzip.open(gtf_path, "rb") as gtf:
    #         for byte_line in gtf:
    #             line = byte_line.decode("utf-8")

    #             if line[0] == "#":
    #                 continue

    #             line = line.rstrip().split("\t")  # type: ignore

    #             # Only process protein lines
    #             if line[2] == "protein":
    #                 # Attempt to fetch both protein_type and protein_biotype, since some ontologies use protein_type and others use protein_biotype
    #                 features = ["protein_id", "protein_name", "protein_version", "protein_type", "protein_biotype"]
    #             else:
    #                 continue

    #             # Extract features (column 9 of GTF)
    #             current_features = gtf_tools._get_features(line)  # type: ignore

    #             # Set protein_name to protein_id if not present in GTF line
    #             if "protein_name" not in current_features:
    #                 current_features["protein_name"] = current_features["protein_id"]

    #             # Filter proteins suffixed with "PAR_Y"
    #             if current_features["protein_id"].endswith("PAR_Y"):
    #                 continue

    #             # get protein length
    #             current_length = protein_lengths[current_features["protein_id"]]

    #             # Select features of interest, raise error if feature of interest not found
    #             target_features = [""] * len(features)
    #             for i in range(len(features)):
    #                 feature = features[i]
    #                 if feature in current_features:
    #                     target_features[i] = current_features[feature]

    #                 # if the symbol starts with ENSG and it does not match the Ensembl ID, then the symbol used should be
    #                 # the Ensembl ID
    #                 if (
    #                     feature in ["protein_name"]
    #                     and current_features[feature].startswith("ENSG")
    #                     and current_features[feature] != current_features["protein_id"]
    #                 ):
    #                     target_features[i] = current_features["protein_id"]

    #                 # Add protein version if available from protein id
    #                 if feature in ["protein_id"]:
    #                     if "." in target_features[i]:
    #                         (feature_id, feature_version) = target_features[i].split(".")
    #                     else:
    #                         feature_id = target_features[i]
    #                         feature_version = ""

    #                     target_features[i] = feature_id
    #                     current_features[feature.replace("id", "version")] = feature_version

    #             protein_id = target_features[0]
    #             self.protein_metadata[protein_id] = ProteinProcessingResult(
    #                 protein_id=target_features[0],
    #                 protein_name=target_features[1],
    #                 protein_version=target_features[2],
    #                 protein_length=str(current_length),
    #                 # Prefer using protein_type, otherwise use protein_biotype
    #                 protein_type=target_features[3] if target_features[3] != "" else target_features[4],
    #             )
    #             if protein_info_description in self.protein_ids_by_description:
    #                 self.protein_ids_by_description[protein_info_description].append(protein_id)
    #             else:
    #                 self.protein_ids_by_description[protein_info_description] = [protein_id]

    # def _get_protein_lengths_from_gtf(self, gtf_path: str) -> Dict[str, int]:
    #     """
    #     Parses a GTF file and calculates protein lengths, which are calculated as follows for each protein:
    #     1. Get lengths for all different isoforms
    #     2. Get the median of the lengths of these isoforms

    #     :param str gtf_path: path to gzipped gtf file

    #     :rtype  Dict[str]
    #     :return A dictionary with keys being protein ids and values the corresponding length in base pairs
    #     """
    #     protein_to_isoforms_map = defaultdict(set)
    #     isoform_to_length_map = defaultdict(int)
    #     with gzip.open(gtf_path, "rb") as gtf:
    #         for byte_line in gtf:
    #             line = byte_line.decode("utf-8")
    #             if line[0] == "#":
    #                 continue

    #             # See https://www.gencodeproteins.org/pages/data_format.html for GTF metadata schema
    #             protein_metadata = line.rstrip().split("\t")  # type: ignore

    #             if protein_metadata[2] != "exon":
    #                 continue

    #             # Calculate exon length using genomic end location and genomic start location
    #             exon_length = int(protein_metadata[4]) - int(protein_metadata[3]) + 1
    #             current_features = gtf_tools._get_features(protein_metadata)  # type: ignore
    #             protein_id = current_features["protein_id"]
    #             transcript_id = current_features["transcript_id"]

    #             protein_to_isoforms_map[protein_id].add(transcript_id)
    #             isoform_to_length_map[transcript_id] += exon_length

    #     protein_lengths = {}
    #     for protein_id in protein_to_isoforms_map:
    #         isoforms = protein_to_isoforms_map[protein_id]
    #         isoform_lengths = []
    #         for isoform in isoforms:
    #             isoform_lengths.append(isoform_to_length_map[isoform])
    #         # GTFTools established standard is to convert to int
    #         protein_lengths[protein_id] = int(statistics.median(isoform_lengths))

    #     return protein_lengths

    # def _process_ercc(self, ercc_path: str, protein_info_description: str):
    #     """
    #     process the ERCC download, keeps only first column with no header

    #     :param str gtf_path: path to ercch download from thermo fisher
    #     :param str output_json_file: path to output file
    #     """
    #     with open(ercc_path, "r") as ercc:
    #         lines = ercc.readlines()[1:]
    #         for line in lines:
    #             line = line.rstrip().split("\t")  # type: ignore
    #             ercc_id = line[0]
    #             errc_length = str(len(line[4]))
    #             errc_version = "1"

    #             self.protein_metadata[ercc_id] = ProteinProcessingResult(
    #                 protein_id=ercc_id,
    #                 protein_name=ercc_id + " (spike-in control)",
    #                 protein_version=errc_version,
    #                 protein_length=errc_length,
    #                 protein_type="synthetic",
    #             )
    #             if protein_info_description in self.protein_ids_by_description:
    #                 self.protein_ids_by_description[protein_info_description].append(ercc_id)
    #             else:
    #                 self.protein_ids_by_description[protein_info_description] = [ercc_id]

    def process_protein_infos(self, protein_infos: dict) -> None:
        for protein_info_key in protein_infos:
            print(protein_info_key)
            # Add to self.protein_labels and self.protein_ids_by_description

            self.process_individual_protein_info(protein_infos[protein_info_key])

        # Deduplicate protein names across all protein_metadata
        # print("Deduplicating protein names...")
        # protein_names = set()
        # duplicated_protein_names = set()
        # for protein_metadata in self.protein_metadata.values():
        #     protein_name = protein_metadata.protein_name
        #     if protein_name in protein_names:
        #         duplicated_protein_names.add(protein_name)
        #     else:
        #         protein_names.add(protein_name)
        # for protein_id, protein_metadata in self.protein_metadata.items():
        #     protein_name = protein_metadata.protein_name
        #     if protein_name in duplicated_protein_names:
        #         self.protein_metadata[protein_id].protein_name = protein_name + "_" + protein_id

        # # Write output for each file, and process file diffs
        # for protein_info_key in protein_infos:
        #     protein_info_description = protein_infos[protein_info_key]["description"]
        #     print("Writing output for", protein_info_description)
        #     new_file = os.path.join(env.GENCODE_DIR, f"new_proteins_{protein_info_description}.csv.gz")
        #     previous_ref_filepath = os.path.join(env.GENCODE_DIR, f"proteins_{protein_info_description}.csv.gz")
        #     protein_ids = self.protein_ids_by_description[protein_info_description]
        #     output_to_print = ""
        #     try:
        #         for protein_id in protein_ids:
        #             protein_metadata = self.protein_metadata[protein_id]
        #             output_to_print += (
        #                 ",".join(
        #                     [
        #                         protein_metadata.protein_id,
        #                         protein_metadata.protein_name,
        #                         protein_metadata.protein_version,
        #                         protein_metadata.protein_length,
        #                         protein_metadata.protein_type,
        #                     ]
        #                 )
        #                 + "\n"
        #             )
        #         self.write_gzip(output_to_print, new_file)

        #         # Rename new_proteins_{filename} to proteins_{filename} if this is the first time we're introducing a new protein file type
        #         if not os.path.exists(previous_ref_filepath):
        #             os.replace(new_file, previous_ref_filepath)

        #         # Process diff between new file and previous file, if there is a difference
        #         else:
        #             if self.digest(new_file) == self.digest(previous_ref_filepath):
        #                 print(
        #                     "New protein reference is identical to previous protein reference", protein_info_description
        #                 )
        #                 os.remove(new_file)
        #             else:
        #                 print("proteinrating protein reference diff for", protein_info_description)
        #                 self.proteinrate_protein_ref_diff(protein_info_description, new_file, previous_ref_filepath)
        #                 os.replace(new_file, previous_ref_filepath)
        #     except Exception as e:
        #         print("Writing to new file failed. Using previous version.", protein_info_description)
        #         raise e

    def process_individual_protein_info(self, protein_info: dict) -> None:
        """
        Download the protein_info and convert it into a csv
        :param protein_info: context to download and process the protein information.
        """
        protein_info_description = protein_info["description"]
        print("Start", protein_info_description)
        url = protein_info["url"]

        print("download", protein_info_description)

        file = os.path.join(env.UNIPROT_DIR, f"proteins_{protein_info_description}.tsv.gz")

        # Download from uniprot without pagination
        with requests.get(url, stream=True) as request:
            request.raise_for_status()
            with open(file, "wb") as f:
                for chunk in request.iter_content(chunk_size=2**20):
                    f.write(chunk)

        # Download from uniprot with pagination
        # re_next_link = re.compile(r'<(.+)>; rel="next"')
        # retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
        # session = requests.Session()
        # session.mount("https://", HTTPAdapter(max_retries=retries))

        # def get_next_link(headers):
        #     if "Link" in headers:
        #         match = re_next_link.match(headers["Link"])
        #         if match:
        #             return match.group(1)

        # def get_batch(batch_url):
        #     while batch_url:
        #         response = session.get(batch_url)
        #         response.raise_for_status()
        #         total = response.headers["x-total-results"]
        #         yield response, total
        #         batch_url = get_next_link(response.headers)

        # progress = 0
        # with open(file, "w") as f:
        #     for batch, total in get_batch(url):
        #         lines = batch.text.splitlines()
        #         if not progress:
        #             print(lines[0], file=f)
        #         for line in lines[1:]:
        #             print(line, file=f)
        #         progress += len(lines[1:])
        #         print(f"{progress} / {total}")
        #         print("finish", protein_info_description)

    # def proteinrate_protein_ref_diff(
    #     self, output_filename: str, current_ref_filepath: str, previous_ref_filepath: str
    # ) -> None:
    #     """
    #     Compare the previous protein reference CSV to the newest proteinrated CSV. Report a text file with every
    #     Ensembl ID available in the previous protein reference CSV that is no longer found in the newest.
    #     :param output_filename: Filename for output diff textfile
    #     :param current_ref_filepath: Full filepath for the new processed protein csv
    #     :param previous_ref_filepath: Full filepath for the previous processed protein csv
    #     """
    #     new_ref_protein_ids = set()
    #     with gzip.open(current_ref_filepath, "rt") as f:
    #         for row in csv.reader(f):
    #             protein_id = row[0]
    #             new_ref_protein_ids.add(protein_id)

    #     removed_protein_ids = []
    #     with gzip.open(previous_ref_filepath, "rt") as f:
    #         for row in csv.reader(f):
    #             protein_id = row[0]
    #             if protein_id not in new_ref_protein_ids:
    #                 removed_protein_ids.append(protein_id)

    #     diff_filepath = os.path.join(env.GENCODE_DIR, f"{output_filename}_diff.txt")
    #     if removed_protein_ids:
    #         with open(diff_filepath, "w") as f:
    #             for protein_id in removed_protein_ids:
    #                 f.write(f"{protein_id}\n")
    #     else:
    #         print(f"No proteins removed from {output_filename} reference. No diff created.")


def main():
    with open(env.PROTEIN_INFO_YAML, "r") as protein_info_handle:
        protein_infos: dict = yaml.safe_load(protein_info_handle)

    protein_processor = ProteinProcessor()
    protein_processor.process_protein_infos(protein_infos)


if __name__ == "__main__":
    main()
