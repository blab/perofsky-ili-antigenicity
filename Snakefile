from augur.frequency_estimators import timestamp_to_float
from datetime import date
import pandas as pd
from treetime.utils import numeric_date

configfile: "config.yaml"
localrules: download_sequences, download_all_titers_by_assay, get_titers_by_passage, filter_metadata

wildcard_constraints:
    sample="sample_(\d+|titers)",
    viruses="\d+",
    bandwidth="[0-9]*\.?[0-9]+",
    lineage="[a-z0-9]+",
    segment="[a-z]+[0-9]?",
    region="[-a-z]+",
    resolution="\d+y",
    replicate="\d+"

path_to_fauna = '../fauna'
segments = ['ha']
lineages = ['h3n2']
resolutions = ['21y']
passages = ['cell']
regions = ["global", "north-america"]
frequency_regions = ['north_america', 'south_america', 'europe', 'china',
                     'southeast_asia', 'japan_korea', 'south_asia', 'africa']
replicates = range(0, config["number_of_replicates"])

def _get_float_date_from_string(date_string):
    return timestamp_to_float(pd.to_datetime(date_string))

MIN_DATE = _get_float_date_from_string(config["start_date"])
MAX_DATE = _get_float_date_from_string(config["end_date"])

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
    return ["results/%s/%s/aa-seq_%s_%s_%s_%s.fasta" % (w.replicate, w.region, w.lineage, w.segment, w.resolution, g)
            for g in genes]

def translations_jsons(w):
    genes = gene_names(w)
    return ["results/%s/%s/aa-seq_%s_%s_%s_%s.json" % (w.replicate, w.region, w.lineage, w.segment, w.resolution, g)
            for g in genes]

def pivots_per_year(w):
    pivots_per_year = {'2y':12, '3y':6, '6y':4, '12y':2}
    return pivots_per_year[w.resolution]

def min_date(w):
    now = numeric_date(date.today())
    return now - int(w.resolution[:-1])

def max_date(w):
    return numeric_date(date.today())

def _get_clock_rate_by_wildcards(wildcards):
    rates_by_lineage_and_segment = {
        ('h3n2', 'ha'): 0.0043, ('h3n2', 'na'):0.0029,
        ('h1n1pdm', 'ha'): 0.0040, ('h1n1pdm', 'na'):0.0032,
        ('vic', 'ha'): 0.0024, ('vic', 'na'):0.0015,
        ('yam', 'ha'): 0.0019, ('yam', 'na'):0.0013
    }

    try:
        rate = rates_by_lineage_and_segment[(wildcards.lineage, wildcards.segment)]
    except KeyError:
        print(f"ERROR: No clock rate defined for {wildcards.lineage} and {wildcards.segment}", file=sys.stderr)
        raise

    return rate

def _get_clock_std_dev_by_wildcards(wildcards):
    return 0.2 * _get_clock_rate_by_wildcards(wildcards)

def vpm(v):
    vpm = {'2y':2, '3y':2, '6y':2, '12y':1, '21y': 50}
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
masks_config = pd.read_table("config/distance_maps.tsv")

def _get_build_mask_config(wildcards):
    config = masks_config[(masks_config["lineage"] == wildcards.lineage) &
                          (masks_config["segment"] == wildcards.segment)]
    if config.shape[0] > 0:
        return config
    else:
        return None

def _get_distance_comparisons_by_lineage_and_segment(wildcards):
    config = _get_build_mask_config(wildcards)
    return " ".join(config.loc[:, "compare_to"].values)

def _get_distance_attributes_by_lineage_and_segment(wildcards):
    config = _get_build_mask_config(wildcards)
    return " ".join(config.loc[:, "attribute"].values)

def _get_seasonal_distance_attributes_by_lineage_and_segment(wildcards):
    config = _get_build_mask_config(wildcards)
    return " ".join(["%s_seasonal" % attribute for attribute in config.loc[:, "attribute"].values])

def _get_distance_maps_by_lineage_and_segment(wildcards):
    config = _get_build_mask_config(wildcards)
    return [
        "config/distance_maps/{wildcards.lineage}/{wildcards.segment}/{distance_map}.json".format(wildcards=wildcards, distance_map=distance_map)
        for distance_map in config.loc[:, "distance_map"].values
    ]

#
# Define rules.
#

rule all:
    input:
        mean_lbi = expand("results/{region}_mean_seasonal_lbi_{lineage}_{segment}_{resolution}.tsv", lineage=lineages, segment=segments, resolution=resolutions, region=regions),
        mean_strain_lbi = expand("results/{region}_strain_seasonal_lbi_{lineage}_{segment}_{resolution}.tsv", lineage=lineages, segment=segments, resolution=resolutions, region=regions),
        mean_distances = expand("results/{region}_mean_seasonal_distances_{lineage}_{segment}_{resolution}.tsv", lineage=lineages, segment=segments, resolution=resolutions, region=regions),
        mean_strain_distances = expand("results/{region}_mean_strain_seasonal_distances_{lineage}_{segment}_{resolution}.tsv", lineage=lineages, segment=segments, resolution=resolutions, region=regions),
        mean_titer_sub_distances = expand("results/{region}_mean_titer_sub_seasonal_distances_{lineage}_{passage}_{segment}_{resolution}.tsv", lineage=lineages, passage=passages, segment=segments, resolution=resolutions, region=regions),
        mean_strain_titer_sub_distances = expand("results/{region}_mean_strain_titer_sub_seasonal_distances_{lineage}_{passage}_{segment}_{resolution}.tsv.gz", lineage=lineages, passage=passages, segment=segments, resolution=resolutions, region=regions),
        mean_titer_tree_distances = expand("results/{region}_mean_titer_tree_seasonal_distances_{lineage}_{passage}_{segment}_{resolution}.tsv", lineage=lineages, passage=passages, segment=segments, resolution=resolutions, region=regions),
        mean_strain_titer_tree_distances = expand("results/{region}_mean_strain_titer_tree_seasonal_distances_{lineage}_{passage}_{segment}_{resolution}.tsv.gz", lineage=lineages, passage=passages, segment=segments, resolution=resolutions, region=regions),
        auspice_tables = expand("auspice_tables/flu_seasonal_{lineage}_{segment}_{resolution}_{region}.tsv", lineage=lineages, segment=segments, resolution=resolutions, region=regions),
        auspice = expand("auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{region}_{replicate}.json", lineage=lineages, segment=segments, resolution=resolutions, region=regions, replicate=replicates),
        auspice_frequencies = expand("auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{region}_{replicate}_tip-frequencies.json", lineage=lineages, segment=segments, resolution=resolutions, region=regions, replicate=replicates)

rule auspice:
    input:
        auspice = expand("auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{region}_{replicate}.json", lineage=lineages, segment=segments, resolution=resolutions, region=regions, replicate=replicates),
        auspice_frequencies = expand("auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{region}_{replicate}_tip-frequencies.json", lineage=lineages, segment=segments, resolution=resolutions, region=regions, replicate=replicates)

rule files:
    params:
        outliers = "config/outliers_{lineage}.txt",
        references = "config/references_{lineage}.txt",
        reference = "config/reference_{lineage}_{segment}.gb",
        colors = "config/colors.tsv",
        auspice_config = "config/auspice_config.json",
        vaccine_json = "config/vaccines_{lineage}.json"

files = rules.files.params

rule download_sequences:
    message: "Downloading sequences from fauna"
    output:
        sequences = "data/{lineage}_{segment}.fasta"
    params:
        fasta_fields = "strain virus accession collection_date region country division location passage_category submitting_lab age gender"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        python3 {path_to_fauna}/vdb/download.py \
                --database vdb \
                --virus flu \
                --fasta_fields {params.fasta_fields} \
                --resolve_method split_passage \
                --select locus:{wildcards.segment} lineage:seasonal_{wildcards.lineage} \
                --path data \
                --fstem {wildcards.lineage}_{wildcards.segment}
        """

rule download_all_titers_by_assay:
    message: "Downloading {wildcards.lineage} HI titers from fauna"
    output:
        titers = "data/{lineage}_hi_titers.json"
    conda: "envs/nextstrain.yaml"
    benchmark: "benchmarks/download_all_titers_{lineage}_hi.txt"
    log: "logs/download_all_titers_{lineage}_hi.log"
    params:
        databases = config["titer_databases"]
    shell:
        """
        python3 {path_to_fauna}/tdb/download.py \
            --database {params.databases} \
            --virus flu \
            --subtype {wildcards.lineage} \
            --select assay_type:hi \
            --path data \
            --fstem {wildcards.lineage}_hi \
            --ftype json
        """

rule get_titers_by_passage:
    message: "Parsing {wildcards.passage}-passaged titers for {wildcards.lineage} HI"
    input:
        titers = rules.download_all_titers_by_assay.output.titers
    output:
        titers = "data/{lineage}_{passage}_hi_titers.tsv"
    benchmark: "benchmarks/get_titers_{lineage}_{passage}_hi.txt"
    log: "logs/get_titers_{lineage}_{passage}_hi.log"
    run:
        df = pd.read_json(input.titers)
        passaged = (df["serum_passage_category"] == wildcards.passage)
        tdb_passaged = df["index"].apply(lambda index: isinstance(index, list) and wildcards.passage in index)
        tsv_fields = [
            "virus_strain",
            "serum_strain",
            "serum_id",
            "source",
            "titer",
            "assay_type"
        ]

        titers_df = df.loc[(passaged | tdb_passaged), tsv_fields]
        titers_df.to_csv(output.titers, sep="\t", header=False, index=False)

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

rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering.
        """
    input:
        sequences = rules.parse.output.sequences
    output:
        sequence_index = "results/sequence_index_{lineage}_{segment}.tsv"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
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
        sequence_index = rules.index_sequences.output.sequence_index,
        include = files.references,
        exclude = files.outliers
    output:
        sequences = 'results/filtered_{lineage}_{segment}.fasta'
    params:
        min_length = config["min_length"]
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --min-length {params.min_length} \
            --non-nucleotide \
            --include {input.include} \
            --exclude {input.exclude} \
            --exclude-where country=? region=? passage=egg \
            --output {output}
        """

rule filter_metadata:
    message:
        """
        Excluding strains with ambiguous dates for {wildcards.lineage} {wildcards.segment}
        """
    input:
        metadata = rules.parse.output.metadata,
        references = files.references
    output:
        metadata = "results/filtered_metadata_{lineage}_{segment}.tsv"
    run:
        references = pd.read_csv(input.references, header=None, names=["strain"])
        df = pd.read_csv(input.metadata, sep="\t")
        is_reference_strain = df["strain"].isin(references["strain"].values)
        has_unambiguous_year_month = ~df["date"].str.contains("XX-XX")
        df[is_reference_strain | has_unambiguous_year_month].to_csv(output.metadata, sep="\t", header=True, index=False)

rule select_strains:
    input:
        sequences = expand("results/filtered_{{lineage}}_{segment}.fasta", segment=segments),
        metadata = expand("results/filtered_metadata_{{lineage}}_{segment}.tsv", segment=segments),
        titers = "data/{lineage}_cell_hi_titers.tsv",
        include = files.references
    output:
        strains = "results/{replicate}/{region}/strains_{lineage}_{resolution}.txt",
    log:
        "logs/strains_{region}_{lineage}_{resolution}_{replicate}.txt"
    params:
        start_date = config["start_date"],
        end_date = config["end_date"],
        viruses_per_month = vpm,
        region = lambda wildcards: wildcards.region.replace("-", "_")
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
            --priority-region {params.region} \
            --output {output.strains} &> {log}
        """

def _get_strains_to_extract(wildcards):
    """Get list of frozen strains for the given wildcards, if they have been
    defined. Otherwise, select strains from scratch.
    """
    frozen_strains = config.get("frozen_strains", {}).get(wildcards.region)

    if frozen_strains:
        return frozen_strains
    else:
        return rules.select_strains.output.strains.format(**wildcards)

rule extract:
    input:
        sequences = rules.filter.output.sequences,
        strains = _get_strains_to_extract
    output:
        sequences = 'results/{replicate}/{region}/extracted_{lineage}_{segment}_{resolution}.fasta'
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
        alignment = "results/{replicate}/{region}/aligned_{lineage}_{segment}_{resolution}.fasta"
    conda: "envs/nextstrain.yaml"
    threads: 8
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --remove-reference \
            --nthreads {threads}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/{replicate}/{region}/tree-raw_{lineage}_{segment}_{resolution}.nwk"
    conda: "envs/nextstrain.yaml"
    threads: 8
    shell:
        """
        export AUGUR_RECURSION_LIMIT=10000 && augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads {threads}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = rules.parse.output.metadata
    output:
        tree = "results/{replicate}/{region}/tree_{lineage}_{segment}_{resolution}.nwk",
        node_data = "results/{replicate}/{region}/branch-lengths_{lineage}_{segment}_{resolution}.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4,
        clock_rate = _get_clock_rate_by_wildcards,
        clock_std_dev = _get_clock_std_dev_by_wildcards
    conda: "envs/nextstrain.yaml"
    log: "logs/refine_{region}_{lineage}_{segment}_{resolution}_{replicate}.txt"
    shell:
        """
        export AUGUR_RECURSION_LIMIT=10000 && augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --no-covariance \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/{replicate}/{region}/nt-muts_{lineage}_{segment}_{resolution}.json"
    params:
        inference = "joint"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/{replicate}/{region}/aa-muts_{lineage}_{segment}_{resolution}.json",
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
        """

rule reconstruct_translations:
    message: "Reconstructing translations required for titer models and frequencies"
    input:
        tree = rules.refine.output.tree,
        node_data = "results/{replicate}/{region}/aa-muts_{lineage}_{segment}_{resolution}.json",
    output:
        aa_alignment = "results/{replicate}/{region}/aa-seq_{lineage}_{segment}_{resolution}_{gene}.fasta"
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
        node_data = "results/{replicate}/{region}/traits_{lineage}_{segment}_{resolution}.json",
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

rule convert_translations_to_json:
    input:
        tree = rules.refine.output.tree,
        translations = "results/{replicate}/{region}/aa-seq_{lineage}_{segment}_{resolution}_{gene}.fasta"
    output:
        translations = "results/{replicate}/{region}/aa-seq_{lineage}_{segment}_{resolution}_{gene}.json"
    shell:
        """
        python3 scripts/convert_translations_to_json.py \
            --tree {input.tree} \
            --alignment {input.translations} \
            --gene-name {wildcards.gene} \
            --attribute-name {wildcards.gene} \
            --output {output.translations}
        """

rule titers_sub:
    input:
        titers = "data/{lineage}_{passage}_hi_titers.tsv",
        aa_muts = rules.translate.output,
        alignments = translations,
        tree = rules.refine.output.tree
    params:
        genes = gene_names
    output:
        titers_model = "results/{replicate}/{region}/titers-sub-model_{lineage}_{passage}_{segment}_{resolution}.json",
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
        titers = "data/{lineage}_{passage}_hi_titers.tsv",
        tree = rules.refine.output.tree
    output:
        titers_model = "results/{replicate}/{region}/titers-tree-model_{lineage}_{passage}_{segment}_{resolution}.json",
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur titers tree \
            --titers {input.titers} \
            --tree {input.tree} \
            --output {output.titers_model}
        """

rule convert_titer_model_to_distance_map:
    input:
        model = rules.titers_sub.output.titers_model
    output:
        distance_map = "results/{replicate}/{region}/titer_substitution_distance_map_{lineage}_{passage}_{segment}_{resolution}.json"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/titer_model_to_distance_map.py \
            --model {input.model} \
            --output {output}
        """

rule clades:
    message: "Annotating clades"
    input:
        tree = rules.refine.output.tree,
        nt_muts = rules.ancestral.output,
        aa_muts = rules.translate.output,
        clade_definitions = "config/clades_{lineage}_{segment}.tsv"
    output:
        clades = "results/{replicate}/{region}/clades_{lineage}_{segment}_{resolution}.json"
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
        distance_maps = _get_distance_maps_by_lineage_and_segment
    params:
        genes = gene_names,
        comparisons = _get_distance_comparisons_by_lineage_and_segment,
        attribute_names = _get_distance_attributes_by_lineage_and_segment
    output:
        distances = "results/{replicate}/{region}/distances_{lineage}_{segment}_{resolution}.json"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur distance \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --compare-to {params.comparisons} \
            --attribute-name {params.attribute_names} \
            --map {input.distance_maps} \
            --output {output}
        """

rule seasonal_distances:
    input:
        tree = rules.refine.output.tree,
        alignments = translations,
        distance_maps = _get_distance_maps_by_lineage_and_segment,
        date_annotations = rules.refine.output.node_data
    params:
        genes = gene_names,
        attribute_names = _get_seasonal_distance_attributes_by_lineage_and_segment,
        start_date = config["distance_start_date"],
        end_date = config["distance_end_date"],
        interval = 12,
        months_prior_to_season = 12
    output:
        distances = "results/{replicate}/{region}/seasonal_distances_{lineage}_{segment}_{resolution}.json"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/seasonal_distance.py \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --attribute-name {params.attribute_names} \
            --map {input.distance_maps} \
            --date-annotations {input.date_annotations} \
            --start-date {params.start_date} \
            --end-date {params.end_date} \
            --interval {params.interval} \
            --months-prior-to-season {params.months_prior_to_season} \
            --output {output}
        """

rule seasonal_titer_distances:
    input:
        tree = rules.refine.output.tree,
        alignments = translations,
        distance_maps = rules.convert_titer_model_to_distance_map.output.distance_map,
        date_annotations = rules.refine.output.node_data
    params:
        genes = gene_names,
        attribute_names = "cTiterSub_seasonal",
        start_date = config["distance_start_date"],
        end_date = config["distance_end_date"],
        interval = 12,
        months_prior_to_season = 12
    output:
        distances = "results/{replicate}/{region}/seasonal_titer_distances_{lineage}_{passage}_{segment}_{resolution}.json"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/seasonal_distance.py \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --attribute-name {params.attribute_names} \
            --map {input.distance_maps} \
            --date-annotations {input.date_annotations} \
            --start-date {params.start_date} \
            --end-date {params.end_date} \
            --interval {params.interval} \
            --months-prior-to-season {params.months_prior_to_season} \
            --output {output}
        """

rule seasonal_titer_tree_distances:
    input:
        tree = rules.refine.output.tree,
        date_annotations = rules.refine.output.node_data,
        titer_tree_model = rules.titers_tree.output.titers_model
    params:
        attribute_names = "cTiter_seasonal",
        start_date = config["distance_start_date"],
        end_date = config["distance_end_date"],
        interval = 12,
        months_prior_to_season = 12
    output:
        distances = "results/{replicate}/{region}/seasonal_titer_tree_distances_{lineage}_{passage}_{segment}_{resolution}.json"
    log:
        "logs/{region}_seasonal_titer_tree_distances_{lineage}_{passage}_{segment}_{resolution}_{replicate}.txt"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/seasonal_distance.py \
            --tree {input.tree} \
            --attribute-name {params.attribute_names} \
            --date-annotations {input.date_annotations} \
            --start-date {params.start_date} \
            --end-date {params.end_date} \
            --interval {params.interval} \
            --months-prior-to-season {params.months_prior_to_season} \
            --titer-tree-model {input.titer_tree_model} \
            --output {output} &> {log}
        """

rule seasonal_vaccine_titer_distances:
    input:
        tree = rules.refine.output.tree,
        alignments = translations,
        distance_maps = rules.convert_titer_model_to_distance_map.output.distance_map,
        date_annotations = rules.refine.output.node_data,
        vaccines = "config/vaccine_start_dates_{vaccine_region}.tsv"
    params:
        genes = gene_names,
        attribute_names = "cTiterSub_seasonal_vaccine",
        start_date = config["distance_start_date"],
        end_date = config["distance_end_date"],
        interval = 12,
        months_prior_to_season = 12
    output:
        distances = "results/{replicate}/{region}/seasonal_vaccine_titer_distances_{lineage}_{passage}_{vaccine_region}_{segment}_{resolution}.json"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/seasonal_distance.py \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --attribute-name {params.attribute_names}_{wildcards.vaccine_region} \
            --map {input.distance_maps} \
            --date-annotations {input.date_annotations} \
            --start-date {params.start_date} \
            --end-date {params.end_date} \
            --interval {params.interval} \
            --months-prior-to-season {params.months_prior_to_season} \
            --vaccines {input.vaccines} \
            --output {output}
        """

rule seasonal_vaccine_titer_tree_distances:
    input:
        tree = rules.refine.output.tree,
        date_annotations = rules.refine.output.node_data,
        titer_tree_model = rules.titers_tree.output.titers_model,
        vaccines = "config/vaccine_start_dates_{vaccine_region}.tsv"
    params:
        attribute_names = "cTiter_seasonal_vaccine",
        start_date = config["distance_start_date"],
        end_date = config["distance_end_date"],
        interval = 12,
        months_prior_to_season = 12
    output:
        distances = "results/{replicate}/{region}/seasonal_vaccine_titer_tree_distances_{lineage}_{passage}_{vaccine_region}_{segment}_{resolution}.json"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/seasonal_distance.py \
            --tree {input.tree} \
            --attribute-name {params.attribute_names}_{wildcards.vaccine_region} \
            --date-annotations {input.date_annotations} \
            --start-date {params.start_date} \
            --end-date {params.end_date} \
            --interval {params.interval} \
            --months-prior-to-season {params.months_prior_to_season} \
            --titer-tree-model {input.titer_tree_model} \
            --vaccines {input.vaccines} \
            --output {output}
        """

rule mean_seasonal_distances:
    input:
        tree = rules.refine.output.tree,
        alignments = translations,
        distance_maps = _get_distance_maps_by_lineage_and_segment,
        date_annotations = rules.refine.output.node_data
    params:
        genes = gene_names,
        attribute_names = _get_seasonal_distance_attributes_by_lineage_and_segment,
        start_date = config["distance_start_date"],
        end_date = config["distance_end_date"],
        interval = 12,
        seasons_away_to_compare = 2
    output:
        distances = "results/{replicate}/{region}_mean_seasonal_distances_{lineage}_{segment}_{resolution}.tsv",
        pairwise_distances = "results/{replicate}/{region}_mean_strain_seasonal_distances_{lineage}_{segment}_{resolution}.tsv"
    log: "logs/mean_seasonal_distances_{region}_{lineage}_{segment}_{resolution}_{replicate}.txt"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/mean_seasonal_distances.py \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --attribute-name {params.attribute_names} \
            --map {input.distance_maps} \
            --date-annotations {input.date_annotations} \
            --start-date {params.start_date} \
            --end-date {params.end_date} \
            --interval {params.interval} \
            --seasons-away-to-compare {params.seasons_away_to_compare} \
            --output {output.distances} \
            --strain-output {output.pairwise_distances} &> {log}
        """

rule mean_titer_sub_seasonal_distances:
    input:
        tree = rules.refine.output.tree,
        alignments = translations,
        distance_maps = rules.convert_titer_model_to_distance_map.output.distance_map,
        date_annotations = rules.refine.output.node_data
    params:
        genes = gene_names,
        attribute_names = "cTiterSub",
        start_date = config["distance_start_date"],
        end_date = config["distance_end_date"],
        interval = 12,
        seasons_away_to_compare = 2
    output:
        distances = "results/{replicate}/{region}_mean_titer_sub_seasonal_distances_{lineage}_{passage}_{segment}_{resolution}.tsv",
        pairwise_distances = "results/{replicate}/{region}_mean_strain_titer_sub_seasonal_distances_{lineage}_{passage}_{segment}_{resolution}.tsv"
    log: "logs/mean_titer_sub_seasonal_distances_{region}_{lineage}_{passage}_{segment}_{resolution}_{replicate}.txt"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/mean_seasonal_distances.py \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --attribute-name {params.attribute_names} \
            --map {input.distance_maps} \
            --date-annotations {input.date_annotations} \
            --start-date {params.start_date} \
            --end-date {params.end_date} \
            --interval {params.interval} \
            --seasons-away-to-compare {params.seasons_away_to_compare} \
            --output {output.distances} \
            --strain-output {output.pairwise_distances} &> {log}
        """


rule mean_titer_tree_seasonal_distances:
    input:
        tree = rules.refine.output.tree,
        date_annotations = rules.refine.output.node_data,
        titer_tree_model = rules.titers_tree.output.titers_model,
    params:
        start_date = config["distance_start_date"],
        end_date = config["distance_end_date"],
        interval = 12,
        seasons_away_to_compare = 2
    output:
        distances = "results/{replicate}/{region}_mean_titer_tree_seasonal_distances_{lineage}_{passage}_{segment}_{resolution}.tsv",
        pairwise_distances = "results/{replicate}/{region}_mean_strain_titer_tree_seasonal_distances_{lineage}_{passage}_{segment}_{resolution}.tsv"
    log: "logs/mean_titer_tree_seasonal_distances_{region}_{lineage}_{passage}_{segment}_{resolution}_{replicate}.txt"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/mean_titer_tree_seasonal_distances.py \
            --tree {input.tree} \
            --date-annotations {input.date_annotations} \
            --model {input.titer_tree_model} \
            --start-date {params.start_date} \
            --end-date {params.end_date} \
            --interval {params.interval} \
            --seasons-away-to-compare {params.seasons_away_to_compare} \
            --output {output.distances} \
            --strain-output {output.pairwise_distances} &> {log}
        """


rule mean_seasonal_lbi:
    input:
        tree = rules.refine.output.tree,
        date_annotations = rules.refine.output.node_data
    params:
        start_date = config["distance_start_date"],
        end_date = config["distance_end_date"],
        interval = 12
    output:
        mean_lbi = "results/{replicate}/{region}_mean_seasonal_lbi_{lineage}_{segment}_{resolution}.tsv",
        lbi_by_strain = "results/{replicate}/{region}_strain_seasonal_lbi_{lineage}_{segment}_{resolution}.tsv"
    log: "logs/mean_seasonal_lbi_{region}_{lineage}_{segment}_{resolution}_{replicate}.txt"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/mean_seasonal_lbi.py \
            --tree {input.tree} \
            --date-annotations {input.date_annotations} \
            --start-date {params.start_date} \
            --end-date {params.end_date} \
            --interval {params.interval} \
            --output {output.mean_lbi} \
            --strain-output {output.lbi_by_strain} &> {log}
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
        lbi = "results/{replicate}/{region}/lbi_{lineage}_{segment}_{resolution}.json"
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

rule tip_frequencies:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        weights = "config/frequency_weights_by_region.json"
    params:
        narrow_bandwidth = 1 / 12.0,
        wide_bandwidth = 3 / 12.0,
        proportion_wide = 0.0,
        weight_attribute = "region",
        min_date = MIN_DATE,
        max_date = MAX_DATE,
        pivot_interval = config["pivot_interval"]
    output:
        tip_freq = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{region}_{replicate}_tip-frequencies.json"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        augur frequencies \
            --method kde \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --narrow-bandwidth {params.narrow_bandwidth} \
            --wide-bandwidth {params.wide_bandwidth} \
            --proportion-wide {params.proportion_wide} \
            --weights {input.weights} \
            --weights-attribute {params.weight_attribute} \
            --pivot-interval {params.pivot_interval} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --output {output}
        """

def _get_node_data_for_export(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by all builds.
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data
    ] + translations_jsons(wildcards)

    # Only run titer rules for segments with titer data.
    if wildcards.segment == "ha":
        inputs.extend([
            rules.titers_tree.output.titers_model,
            rules.titers_sub.output.titers_model
        ])

    # Only request a distance file for builds that have mask configurations
    # defined.
    if _get_build_mask_config(wildcards) is not None:
        inputs.append(rules.distances.output.distances)
        inputs.append(rules.seasonal_distances.output.distances)

        if wildcards.segment == "ha":
            inputs.append(rules.seasonal_titer_distances.output.distances)
            inputs.append(rules.seasonal_titer_tree_distances.output.distances)
            inputs.append(rules.seasonal_vaccine_titer_distances.output.distances)
            inputs.append(rules.seasonal_vaccine_titer_tree_distances.output.distances)

    # Convert input files from wildcard strings to real file names.
    wildcards_dict = dict(wildcards)
    wildcards_dict["passage"] = "cell"
    wildcards_dict["vaccine_region"] = "northern-hemisphere"

    inputs = [input_file.format(**wildcards_dict) for input_file in inputs]
    return inputs

rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        auspice_config = files.auspice_config,
        vaccines = files.vaccine_json,
        node_data = _get_node_data_for_export
    output:
        auspice_main = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{region}_{replicate}.json"
    params:
        color_by_fields = "date"
    conda: "envs/nextstrain.yaml"
    shell:
        """
        export AUGUR_RECURSION_LIMIT=10000 && augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.vaccines} {input.node_data} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_main} \
            --minify-json \
            --include-root-sequence \
            --color-by-metadata {params.color_by_fields}
        """


def _get_attributes_to_export(wildcards):
    return config["attributes_to_export_%s" % wildcards.segment]


rule convert_tree_to_table:
    input:
        tree = rules.export.output.auspice_main
    output:
        table = "auspice_tables/{replicate}/flu_seasonal_{lineage}_{segment}_{resolution}_{region}.tsv"
    params:
        attributes = _get_attributes_to_export
    conda: "envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/tree_to_table.py \
            {input.tree} \
            {output} \
            --attributes {params.attributes} \
            --annotations replicate={wildcards.replicate}
        """


rule aggregate_mean_seasonal_lbi:
    input:
        tables=expand("results/{replicate}/{{region}}_mean_seasonal_lbi_{{lineage}}_{{segment}}_{{resolution}}.tsv", replicate=replicates),
    output:
        table="results/{region}_mean_seasonal_lbi_{lineage}_{segment}_{resolution}.tsv",
    conda: "envs/nextstrain.yaml"
    params:
        replicates=list(replicates),
    shell:
        """
        python3 scripts/concatenate_tables.py \
            --tables {input.tables} \
            --ids {params.replicates} \
            --output {output.table}
        """


rule aggregate_strain_seasonal_lbi:
    input:
        tables=expand("results/{replicate}/{{region}}_strain_seasonal_lbi_{{lineage}}_{{segment}}_{{resolution}}.tsv", replicate=replicates),
    output:
        table="results/{region}_strain_seasonal_lbi_{lineage}_{segment}_{resolution}.tsv",
    conda: "envs/nextstrain.yaml"
    params:
        replicates=list(replicates),
    shell:
        """
        python3 scripts/concatenate_tables.py \
            --tables {input.tables} \
            --ids {params.replicates} \
            --output {output.table}
        """


rule aggregate_mean_seasonal_distances:
    input:
        tables=expand("results/{replicate}/{{region}}_mean_seasonal_distances_{{lineage}}_{{segment}}_{{resolution}}.tsv", replicate=replicates),
    output:
        table="results/{region}_mean_seasonal_distances_{lineage}_{segment}_{resolution}.tsv",
    conda: "envs/nextstrain.yaml"
    params:
        replicates=list(replicates),
    shell:
        """
        python3 scripts/concatenate_tables.py \
            --tables {input.tables} \
            --ids {params.replicates} \
            --output {output.table}
        """


rule aggregate_mean_strain_seasonal_distances:
    input:
        tables=expand("results/{replicate}/{{region}}_mean_strain_seasonal_distances_{{lineage}}_{{segment}}_{{resolution}}.tsv", replicate=replicates),
    output:
        table="results/{region}_mean_strain_seasonal_distances_{lineage}_{segment}_{resolution}.tsv",
    conda: "envs/nextstrain.yaml"
    params:
        replicates=list(replicates),
    shell:
        """
        python3 scripts/concatenate_tables.py \
            --tables {input.tables} \
            --ids {params.replicates} \
            --output {output.table}
        """


rule aggregate_mean_titer_sub_seasonal_distances:
    input:
        tables=expand("results/{replicate}/{{region}}_mean_titer_sub_seasonal_distances_{{lineage}}_{{passage}}_{{segment}}_{{resolution}}.tsv", replicate=replicates),
    output:
        table="results/{region}_mean_titer_sub_seasonal_distances_{lineage}_{passage}_{segment}_{resolution}.tsv",
    conda: "envs/nextstrain.yaml"
    params:
        replicates=list(replicates),
    shell:
        """
        python3 scripts/concatenate_tables.py \
            --tables {input.tables} \
            --ids {params.replicates} \
            --output {output.table}
        """


rule aggregate_mean_strain_titer_sub_seasonal_distances:
    input:
        tables=expand("results/{replicate}/{{region}}_mean_strain_titer_sub_seasonal_distances_{{lineage}}_{{passage}}_{{segment}}_{{resolution}}.tsv", replicate=replicates),
    output:
        table="results/{region}_mean_strain_titer_sub_seasonal_distances_{lineage}_{passage}_{segment}_{resolution}.tsv.gz",
    conda: "envs/nextstrain.yaml"
    params:
        replicates=list(replicates),
    shell:
        """
        python3 scripts/concatenate_tables.py \
            --tables {input.tables} \
            --ids {params.replicates} \
            --output {output.table}
        """


rule aggregate_mean_titer_tree_seasonal_distances:
    input:
        tables=expand("results/{replicate}/{{region}}_mean_titer_tree_seasonal_distances_{{lineage}}_{{passage}}_{{segment}}_{{resolution}}.tsv", replicate=replicates),
    output:
        table="results/{region}_mean_titer_tree_seasonal_distances_{lineage}_{passage}_{segment}_{resolution}.tsv",
    conda: "envs/nextstrain.yaml"
    params:
        replicates=list(replicates),
    shell:
        """
        python3 scripts/concatenate_tables.py \
            --tables {input.tables} \
            --ids {params.replicates} \
            --output {output.table}
        """

rule aggregate_mean_strain_titer_tree_seasonal_distances:
    input:
        tables=expand("results/{replicate}/{{region}}_mean_strain_titer_tree_seasonal_distances_{{lineage}}_{{passage}}_{{segment}}_{{resolution}}.tsv", replicate=replicates),
    output:
        table="results/{region}_mean_strain_titer_tree_seasonal_distances_{lineage}_{passage}_{segment}_{resolution}.tsv.gz",
    conda: "envs/nextstrain.yaml"
    params:
        replicates=list(replicates),
    shell:
        """
        python3 scripts/concatenate_tables.py \
            --tables {input.tables} \
            --ids {params.replicates} \
            --output {output.table}
        """

rule aggregate_auspice_tables:
    input:
        tables=expand("auspice_tables/{replicate}/flu_seasonal_{{lineage}}_{{segment}}_{{resolution}}_{{region}}.tsv", replicate=replicates),
    output:
        table="auspice_tables/flu_seasonal_{lineage}_{segment}_{resolution}_{region}.tsv",
    conda: "envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/concatenate_tables.py \
            --tables {input.tables} \
            --output {output.table}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
