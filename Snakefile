from datetime import date
import pandas as pd
from treetime.utils import numeric_date

configfile: "config.json"
localrules: download_sequences, download_titers

path_to_fauna = '../fauna'
segments = ['ha']
lineages = ['h3n2']
resolutions = ['21y']
frequency_regions = ['north_america', 'south_america', 'europe', 'china',
                     'southeast_asia', 'japan_korea', 'south_asia', 'africa']


def reference_strain(v):
    references = {'h3n2':"A/Beijing/32/1992",
                  'h1n1pdm':"A/California/07/2009",
                  'vic':"B/HongKong/02/1993",
                  'yam':"B/Singapore/11/1994"
                  }
    return references[v.lineage]

genes_to_translate = {'ha':['SigPep', 'HA1', 'HA2'], 'na':['NA']}
def gene_names(w):
    return genes_to_translate[w.segment]

def translations(w):
    genes = gene_names(w)
    return ["results/aa-seq_%s_%s_%s_%s.fasta"%(w.lineage, w.segment, w.resolution, g)
            for g in genes]

def pivots_per_year(w):
    pivots_per_year = {'2y':12, '3y':6, '6y':4, '12y':2}
    return pivots_per_year[w.resolution]

def min_date(w):
    now = numeric_date(date.today())
    return now - int(w.resolution[:-1])

def max_date(w):
    return numeric_date(date.today())

def substitution_rates(w):
    references = {('h3n2', 'ha'): 0.0038, ('h3n2', 'na'):0.0028}
    return references[(w.lineage, w.segment)]

def vpm(v):
    vpm = {'2y':2, '3y':2, '6y':2, '12y':1, '21y': 36}
    return vpm[v.resolution] if v.resolution in vpm else 5

#
# Define LBI parameters and functions.
#
LBI_params = {
    '2y': {"tau": 0.3, "time_window": 0.5},
    '3y': {"tau": 0.4, "time_window": 0.6},
    '6y': {"tau": 0.25, "time_window": 0.75},
    '12y': {"tau": 0.25, "time_window": 0.75}
}

def _get_lbi_tau_for_wildcards(wildcards):
    return LBI_params[wildcards.resolution]["tau"]

def _get_lbi_window_for_wildcards(wildcards):
    return LBI_params[wildcards.resolution]["time_window"]

#
# Configure amino acid distance masks.
#

# Load mask configuration including which masks map to which attributes per
# lineage and segment.
masks_config = pd.read_table("config/mask_config.tsv")

def _get_build_mask_config(wildcards):
    config = masks_config[(masks_config["lineage"] == wildcards.lineage) &
                          (masks_config["segment"] == wildcards.segment)]
    if config.shape[0] > 0:
        return config
    else:
        return None

def _get_mask_attribute_names_by_wildcards(wildcards):
    config = _get_build_mask_config(wildcards)
    return " ".join(config.loc[:, "attribute"].values)

def _get_mask_names_by_wildcards(wildcards):
    config = _get_build_mask_config(wildcards)
    return " ".join(config.loc[:, "mask"].values)

#
# Define rules.
#

rule all:
    input:
        auspice_tables = expand("auspice_tables/flu_seasonal_{lineage}_{segment}_{resolution}.tsv", lineage=lineages, segment=segments, resolution=resolutions),
        auspice = expand("auspice/flu_seasonal_{lineage}_{segment}_{resolution}.json", lineage=lineages, segment=segments, resolution=resolutions)

rule files:
    params:
        outliers = "config/outliers_{lineage}.txt",
        references = "config/references_{lineage}.txt",
        reference = "config/reference_{lineage}_{segment}.gb",
        colors = "config/colors.tsv",
        auspice_config = "config/auspice_config.json",

files = rules.files.params

rule download_sequences:
    message: "Downloading sequences from fauna"
    output:
        sequences = "data/{lineage}_{segment}.fasta"
    params:
        fasta_fields = "strain virus accession collection_date region country division location passage_category submitting_lab age gender"
    conda: "envs/nextstrain.python2.yaml"
    shell:
        """
        env PYTHONPATH={path_to_fauna} \
            python2 {path_to_fauna}/vdb/download.py \
                --database vdb \
                --virus flu \
                --fasta_fields {params.fasta_fields} \
                --resolve_method split_passage \
                --select locus:{wildcards.segment} lineage:seasonal_{wildcards.lineage} \
                --path data \
                --fstem {wildcards.lineage}_{wildcards.segment}
        """

rule download_titers:
    message: "Downloading titers from fauna"
    output:
        titers = "data/{lineage}_hi_titers.tsv"
    params:
        fasta_fields = "strain virus accession collection_date region country division location passage_category submitting_lab age gender"
    conda: "envs/nextstrain.python2.yaml"
    shell:
        """
        env PYTHONPATH={path_to_fauna} \
            python2 {path_to_fauna}/tdb/download.py \
                --database tdb cdc_tdb \
                --virus flu \
                --subtype {wildcards.lineage} \
                --select assay_type:hi \
                --path data \
                --fstem {wildcards.lineage}_hi
        """

rule parse:
    message: "Parsing fasta into sequences and metadata"
    input:
        sequences = rules.download_sequences.output.sequences
    output:
        sequences = "results/sequences_{lineage}_{segment}.fasta",
        metadata = "results/metadata_{lineage}_{segment}.tsv"
    params:
        fasta_fields =  "strain virus isolate_id date region country division location passage authors age gender"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields}
        """

rule filter:
    message:
        """
        Filtering {wildcards.lineage} {wildcards.segment} sequences:
          - less than {params.min_length} bases
          - outliers
          - samples with missing region and country metadata
        """
    input:
        metadata = rules.parse.output.metadata,
        sequences = rules.parse.output.sequences,
        exclude = files.outliers
    output:
        sequences = 'results/filtered_{lineage}_{segment}.fasta'
    params:
        min_length = config["min_length"]
    conda: "envs/nextstrain.yaml"
    shell:
        # TODO: Exclude egg strains and require vaccine strains
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --min-length {params.min_length} \
            --non-nucleotide \
            --exclude {input.exclude} \
            --exclude-where country=? region=? \
            --output {output}
        """

rule select_strains:
    input:
        sequences = expand("results/filtered_{{lineage}}_{segment}.fasta", segment=segments),
        metadata = expand("results/metadata_{{lineage}}_{segment}.tsv", segment=segments),
        titers = rules.download_titers.output.titers,
        include = files.references
    output:
        strains = "results/strains_{lineage}_{resolution}.txt",
    params:
        start_date = config["start_date"],
        end_date = config["end_date"],
        viruses_per_month = vpm
    conda: "envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/select_strains.py \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --segments {segments} \
            --include {input.include} \
            --lineage {wildcards.lineage} \
            --time-interval {params.start_date} {params.end_date} \
            --viruses_per_month {params.viruses_per_month} \
            --titers {input.titers} \
            --output {output.strains}
        """

rule extract:
    input:
        sequences = rules.filter.output.sequences,
        strains = rules.select_strains.output.strains
    output:
        sequences = 'results/extracted_{lineage}_{segment}_{resolution}.fasta'
    conda: "envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/extract_sequences.py \
            --sequences {input.sequences} \
            --samples {input.strains} \
            --output {output}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.extract.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned_{lineage}_{segment}_{resolution}.fasta"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --remove-reference
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree-raw_{lineage}_{segment}_{resolution}.nwk"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = rules.parse.output.metadata
    output:
        tree = "results/tree_{lineage}_{segment}_{resolution}.nwk",
        node_data = "results/branch-lengths_{lineage}_{segment}_{resolution}.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt-muts_{lineage}_{segment}_{resolution}.json"
    params:
        inference = "joint"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/aa-muts_{lineage}_{segment}_{resolution}.json",
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """

rule reconstruct_translations:
    message: "Reconstructing translations required for titer models and frequencies"
    input:
        tree = rules.refine.output.tree,
        node_data = "results/aa-muts_{lineage}_{segment}_{resolution}.json",
    output:
        aa_alignment = "results/aa-seq_{lineage}_{segment}_{resolution}_{gene}.fasta"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur reconstruct-sequences \
            --tree {input.tree} \
            --mutations {input.node_data} \
            --gene {wildcards.gene} \
            --output {output.aa_alignment} \
            --internal-nodes
        """

rule traits:
    message:
        """
        Inferring ancestral traits for {params.columns!s}
        """
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata
    output:
        node_data = "results/traits_{lineage}_{segment}_{resolution}.json",
    params:
        columns = "region"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

rule titers_sub:
    input:
        titers = rules.download_titers.output.titers,
        aa_muts = rules.translate.output,
        alignments = translations,
        tree = rules.refine.output.tree
    params:
        genes = gene_names
    output:
        titers_model = "results/titers-sub-model_{lineage}_{segment}_{resolution}.json",
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur titers sub \
            --titers {input.titers} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --tree {input.tree} \
            --output {output.titers_model}
        """

rule titers_tree:
    input:
        titers = rules.download_titers.output.titers,
        tree = rules.refine.output.tree
    output:
        titers_model = "results/titers-tree-model_{lineage}_{segment}_{resolution}.json",
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur titers tree \
            --titers {input.titers} \
            --tree {input.tree} \
            --output {output.titers_model}
        """

rule mutation_frequencies:
    input:
        metadata = rules.parse.output.metadata,
        alignment = translations
    params:
        genes = gene_names,
        min_date = min_date,
        max_date = max_date,
        pivots_per_year = pivots_per_year
    output:
        mut_freq = "results/mutation-frequencies_{lineage}_{segment}_{resolution}.json"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur frequencies \
            --alignments {input.alignment} \
            --metadata {input.metadata} \
            --gene-names {params.genes} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --pivots-per-year {params.pivots_per_year} \
            --output {output.mut_freq}
        """

rule tree_frequencies:
    input:
        metadata = rules.parse.output.metadata,
        tree = rules.refine.output.tree
    params:
        regions = frequency_regions + ['global'],
        min_date = min_date,
        max_date = max_date,
        pivots_per_year = pivots_per_year
    output:
        tree_freq = "results/tree-frequencies_{lineage}_{segment}_{resolution}.json",
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur frequencies \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --pivots-per-year {params.pivots_per_year} \
            --regions {params.regions} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --output {output.tree_freq}
        """

rule clades:
    message: "Annotating clades"
    input:
        tree = rules.refine.output.tree,
        nt_muts = rules.ancestral.output,
        aa_muts = rules.translate.output,
        clade_definitions = "config/clades_{lineage}_{segment}.tsv"
    output:
        clades = "results/clades_{lineage}_{segment}_{resolution}.json"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.nt_muts} {input.aa_muts} \
            --clades {input.clade_definitions} \
            --output {output.clades}
        """

rule distances:
    input:
        tree = rules.refine.output.tree,
        alignments = translations,
        masks = "config/{segment}_masks.tsv"
    params:
        genes = gene_names,
        attribute_names = _get_mask_attribute_names_by_wildcards,
        mask_names = _get_mask_names_by_wildcards
    output:
        distances = "results/distances_{lineage}_{segment}_{resolution}.json",
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur distance \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --masks {input.masks} \
            --output {output} \
            --attribute-names {params.attribute_names} \
            --mask-names {params.mask_names}
        """

rule lbi:
    message: "Calculating LBI"
    input:
        tree = rules.refine.output.tree,
        branch_lengths = rules.refine.output.node_data
    params:
        tau = _get_lbi_tau_for_wildcards,
        window = _get_lbi_window_for_wildcards,
        names = "lbi"
    output:
        lbi = "results/lbi_{lineage}_{segment}_{resolution}.json"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur lbi \
            --tree {input.tree} \
            --branch-lengths {input.branch_lengths} \
            --output {output} \
            --attribute-names {params.names} \
            --tau {params.tau} \
            --window {params.window}
        """

def _get_node_data_for_export(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by all builds.
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.titers_tree.output.titers_model,
        rules.titers_sub.output.titers_model
    ]

    # Only request a distance file for builds that have mask configurations
    # defined.
    if _get_build_mask_config(wildcards) is not None:
        inputs.append(rules.distances.output.distances)

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards) for input_file in inputs]
    return inputs

rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        auspice_config = files.auspice_config,
        node_data = _get_node_data_for_export
    output:
        auspice_tree = "auspice_split/flu_seasonal_{lineage}_{segment}_{resolution}_tree.json",
        auspice_meta = "auspice_split/flu_seasonal_{lineage}_{segment}_{resolution}_meta.json"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_meta}
        """

rule convert_tree_to_table:
    input:
        tree = rules.export.output.auspice_tree
    output:
        table = "auspice_tables/flu_seasonal_{lineage}_{segment}_{resolution}.tsv"
    params:
        attributes = config["attributes_to_export"]
    conda: "envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/tree_to_table.py \
            {input.tree} \
            {output} \
            --attributes {params.attributes}
        """

rule merge_auspice_jsons:
    input:
        tree = rules.export.output.auspice_tree,
        metadata = rules.export.output.auspice_meta
    output:
        table = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}.json"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/merge_auspice_jsons.py \
            {input.tree} \
            {input.metadata} \
            {output}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
