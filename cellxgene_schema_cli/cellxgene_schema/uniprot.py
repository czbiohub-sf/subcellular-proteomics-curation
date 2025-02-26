import enum
import gzip
import os
from typing import Union

import pandas as pd

from . import env


class SupportedOrganisms(enum.Enum):
    # NOTE: these could be enumerated from loading the `schema_definition.yaml` and scraping the 'organism_ontology_term_id' constraints
    HOMO_SAPIENS = "NCBITaxon:9606"
    # MUS_MUSCULUS = "NCBITaxon:10090"
    # SARS_COV_2 = "NCBITaxon:2697049"
    # ERCC = "NCBITaxon:32630"
    # DROSOPHILA_MELANOGASTER = "NCBITaxon:7227"
    # DANIO_RERIO = "NCBITaxon:7955"
    # CAENORHABDITIS_ELEGANS = "NCBITaxon:6239"
    # MACACA_FASCICULARIS = "NCBITaxon:9541"
    # ORYCTOLAGUS_CUNICULUS = "NCBITaxon:9986"
    # CALLITHRIX_JACCHUS = "NCBITaxon:9483"
    # GORILLA_GORILLA = "NCBITaxon:9595"
    # MACACA_MULATTA = "NCBITaxon:9544"
    # PAN_TROGLODYTES = "NCBITaxon:9598"
    # SUS_SCROFA = "NCBITaxon:9823"
    # MICROCEBUS_MURINUS = "NCBITaxon:30608"
    # RATTUS_NORVEGICUS = "NCBITaxon:10116"


def get_organism_from_feature_id(
    feature_id: str,
) -> Union[SupportedOrganisms, None]:
    """
    Determines organism based on which gene file the feature id was in

    :param str feature_id: the feature id

    :rtype Union[ontology.SypportedOrganisms, None]
    :return: the organism the feature id is from
    """

    for organism in SupportedOrganisms:
        gene_checker = get_gene_checker(organism)
        if gene_checker.is_valid_id(feature_id):
            return organism

    return None


class GeneChecker:
    """Handles checking gene ids, retrieves symbols"""

    GENE_FILES = {
        SupportedOrganisms.HOMO_SAPIENS: os.path.join(env.UNIPROT_DIR, "proteins_homo_sapiens.tsv.gz"),
        # SupportedOrganisms.MUS_MUSCULUS: os.path.join(env.UNIPROT_DIR, "proteins_mus_musculus.tsv.gz"),
        # SupportedOrganisms.SARS_COV_2: os.path.join(env.UNIPROT_DIR, "proteins_sars_cov_2.tsv.gz"),
        # SupportedOrganisms.ERCC: os.path.join(env.UNIPROT_DIR, "proteins_ercc.tsv.gz"),
        # SupportedOrganisms.DROSOPHILA_MELANOGASTER: os.path.join(
        #     env.UNIPROT_DIR, "proteins_drosophila_melanogaster.tsv.gz"
        # ),
        # SupportedOrganisms.DANIO_RERIO: os.path.join(env.UNIPROT_DIR, "proteins_danio_rerio.tsv.gz"),
        # SupportedOrganisms.CAENORHABDITIS_ELEGANS: os.path.join(env.UNIPROT_DIR, "proteins_caenorhabditis_elegans.tsv.gz"),
        # SupportedOrganisms.MACACA_FASCICULARIS: os.path.join(env.GENCODE_DIR, "genes_macaca_fascicularis.csv.gz"),
        # SupportedOrganisms.ORYCTOLAGUS_CUNICULUS: os.path.join(env.GENCODE_DIR, "genes_oryctolagus_cuniculus.csv.gz"),
        # SupportedOrganisms.CALLITHRIX_JACCHUS: os.path.join(env.GENCODE_DIR, "genes_callithrix_jacchus.csv.gz"),
        # SupportedOrganisms.GORILLA_GORILLA: os.path.join(env.GENCODE_DIR, "genes_gorilla_gorilla.csv.gz"),
        # SupportedOrganisms.MACACA_MULATTA: os.path.join(env.GENCODE_DIR, "genes_macaca_mulatta.csv.gz"),
        # SupportedOrganisms.PAN_TROGLODYTES: os.path.join(env.GENCODE_DIR, "genes_pan_troglodytes.csv.gz"),
        # SupportedOrganisms.SUS_SCROFA: os.path.join(env.GENCODE_DIR, "genes_sus_scrofa.csv.gz"),
        # SupportedOrganisms.MICROCEBUS_MURINUS: os.path.join(env.GENCODE_DIR, "genes_microcebus_murinus.csv.gz"),
        # SupportedOrganisms.RATTUS_NORVEGICUS: os.path.join(env.GENCODE_DIR, "genes_rattus_norvegicus.csv.gz"),
    }

    def __init__(self, species: SupportedOrganisms):
        """
        :param enum.Enum.SupportedSpecies species: item from SupportedOrganisms
        """

        if species not in self.GENE_FILES:
            raise ValueError(f"{species} not supported.")

        self.species = species
        self.gene_dict = {}
        uniprot_df = pd.read_csv(self.GENE_FILES[species], sep="\t")
        for index, row in uniprot_df.iterrows():
            gene_id = row["Entry"]
            gene_label = row["Gene Names (primary)"]
            gene_length = row["Length"]
            gene_location = row["Subcellular location [CC]"]
            self.gene_dict[gene_id] = (gene_label, gene_length, gene_location)

        # with gzip.open(self.GENE_FILES[species], "rt") as genes:
        #     for gene in genes:
        #         gene = gene.rstrip().split("\t")  # type: ignore
        #         if gene[0] == "Entry":
        #             continue
        #         gene_id = gene[0]
        #         gene_label = gene[1]
        #         gene_length = int(gene[3])
        #         gene_type = gene[4]

        #         self.gene_dict[gene_id] = (gene_label, gene_length, gene_type)

    def is_valid_id(self, gene_id: str) -> bool:
        """
        Checks for validity of gene id

        :param str gene_id: ENSEMBL gene id

        :rtype bool
        :return True if the gene_id is a valid ENSEMBL id, False otherwise
        """

        return gene_id in self.gene_dict

    def get_symbol(self, gene_id: str) -> str:
        """
        Gets symbol associated to the ENSEBML id

        :param str gene_id: ENSEMBL gene id

        :rtype str
        :return A gene symbol
        """

        if self.is_valid_id(gene_id):
            return self.gene_dict[gene_id][0]
        else:
            raise ValueError(f"The id '{gene_id}' is not a valid ENSEMBL id for '{self.species}'")

    def get_length(self, gene_id: str) -> int:
        """
        Gets feature length associated to the ENSEBML id

        :param str gene_id: ENSEMBL gene id

        :rtype int
        :return A gene length
        """

        if self.is_valid_id(gene_id):
            return self.gene_dict[gene_id][1]
        else:
            raise ValueError(f"The id '{gene_id}' is not a valid ENSEMBL id for '{self.species}'")

    def get_location(self, gene_id: str) -> str:
        """
        Gets feature type associated to the ENSEBML id

        :param str gene_id: ENSEMBL gene id

        :rtype str
        :return A feature type
        """

        if self.is_valid_id(gene_id):
            return self.gene_dict[gene_id][2]
        else:
            raise ValueError(f"The id '{gene_id}' is not a valid ENSEMBL id for '{self.species}'")


# cache the gene checkers
_gene_checkers = {}


def get_gene_checker(species: SupportedOrganisms) -> GeneChecker:
    # Values will be instances of gencode.GeneChecker,
    # keys will be one of gencode.SupportedOrganisms
    if species not in _gene_checkers:
        _gene_checkers[species] = GeneChecker(species)
    return _gene_checkers[species]
