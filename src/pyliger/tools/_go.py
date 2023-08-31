from pathlib import Path

import mygene
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.base import download_go_basic_obo, download_ncbi_associations
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.obo_parser import GODag

# from __future__ import print_function


def run_GO_analysis(
    gene_list, background, data_source, result_path=None, alpha=0.05, methods=["fdr_bh"]
):
    """
    Wrapper function to run GOATOOLS
    :param gene_list:
    :param background:
    :param data_source:
    alpha: significance cut-off, default 0.5
    methods: defult multipletest correction method
    :return:
    """
    ### 0. Preprocessing parameters
    # determine data source
    if data_source == "human":
        taxids = [9606]
    elif data_source == "mouse":
        taxids = [10090]
    else:
        print("Data source can only be either human or mouse.")

    if result_path is None:
        obo = Path("./results", "go-basic.obo")
        gene2go = Path("./results", "gene2go")
    else:
        obo = Path(result_path, "go-basic.obo")
        gene2go = Path(result_path, "gene2go")

    ### 1. Download Ontologies and Associations
    obo_fname, fin_gene2go = _download_go(obo, gene2go)

    ### 2a. Load Ontologies and Associations
    obodag, ns2assoc = _load_go(obo_fname, fin_gene2go, taxids)

    ### 2b. Process background gene set
    if type(background[0]) == int:
        background_ids = background
    elif type(background[0]) == str:
        background_ids = _symbols_to_ids(background, taxids[0])

    ### 3. Initialize a GOEA object
    goeaobj = GOEnrichmentStudyNS(
        background_ids,  # List of mouse protein-coding genes
        ns2assoc,  # geneid/GO associations
        obodag,  # Ontologies
        propagate_counts=False,
        alpha=alpha,  # default significance cut-off
        methods=methods,
    )  # defult multipletest correction method

    ### 4. Load study genes
    if type(gene_list[0]) == int:
        geneids_study = gene_list
    elif type(gene_list[0]) == str:
        geneids_study = _symbols_to_ids(gene_list, taxids[0])

    ### 5. Run GOEA
    goea_results_all = goeaobj.run_study(geneids_study)
    # goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]

    return goea_results_all


def _symbols_to_ids(symbols, species):
    """Convert gene symbols to gene ids using mygene module"""
    mg = mygene.MyGeneInfo()
    mg_df = mg.querymany(symbols, scopes="symbol", species=species, as_dataframe=True)
    ids = mg_df["_id"].dropna()
    ids = ids[~ids.str.startswith("ENSG")]
    ids = ids.astype(int)
    ids = list(ids)

    return ids


def _download_go(obo, gene2go):
    """helper function to download Ontologies and Associations"""
    # Ontologies
    obo_fname = download_go_basic_obo(obo=obo)
    # Associations
    fin_gene2go = download_ncbi_associations(gene2go=gene2go)

    return obo_fname, fin_gene2go


def _load_go(obo_fname, fin_gene2go, taxids):
    """helper function to load Ontologies and Associations"""
    # Ontologies
    obodag = GODag(obo_fname)

    # Associaions
    # Read NCBI's gene2go. Store annotations in a list of namedtuples
    objanno = Gene2GoReader(fin_gene2go, taxids=taxids)

    # Get namespace2association where:
    #    namespace is:
    #        BP: biological_process
    #        MF: molecular_function
    #        CC: cellular_component
    #    assocation is a dict:
    #        key: NCBI GeneID
    #        value: A set of GO IDs associated with that gene
    ns2assoc = objanno.get_ns2assc()

    return obodag, ns2assoc
