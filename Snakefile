path_to_fauna = '../fauna/data/seasonal-flu'
segments = ['ha']
lineages = ['h3n2']
resolutions = ['2y']

def reference_strain(v):
    references = {'h3n2':"A/Beijing/32/1992",
                  'h1n1pdm':"A/California/07/2009",
                  'vic':"B/HongKong/02/1993",
                  'yam':"B/Singapore/11/1994"
                  }
    return references[v.lineage]

def titer_data(w):
    titers = {'h1n1':path_to_fauna + '/h1n1_crick_hi_cell_titers.tsv',
            'h3n2':path_to_fauna + '/h3n2_crick_hi_cell_titers.tsv',
            'yam':path_to_fauna + '/yam_crick_hi_cell_titers.tsv',
            'vic':path_to_fauna + '/vic_crick_hi_cell_titers.tsv',
            'Ball':path_to_fauna + '/Ball_crick_hi_cell_titers.tsv',
            'h1n1pdm':path_to_fauna + '/h1n1pdm_crick_hi_cell_titers.tsv'}
    return titers[w.lineage]

def priority_files(w):
    priorty = {'h1n1':path_to_fauna + '/h1n1_crick_hi_cell_strains.tsv',
            'h3n2':path_to_fauna + '/h3n2_crick_hi_cell_strains.tsv',
            'yam':path_to_fauna + '/yam_crick_hi_cell_strains.tsv',
            'vic':path_to_fauna + '/vic_crick_hi_cell_strains.tsv',
            'Ball':path_to_fauna + '/Ball_crick_hi_cell_strains.tsv',
            'h1n1pdm':path_to_fauna + '/h1n1pdm_crick_hi_cell_strains.tsv'}

    return priorty[w.lineage]


def vpm(v):
    vpm = {'3y':2, '6y':2, '12y':1}
    return vpm[v.resolution] if v.resolution in vpm else 5


rule all:
    input:
        auspice_tree = "auspice/flu_{build}_{lineage}_{segment}_{resolution}_tree.json",
        auspice_meta = "auspice/flu_{build}_{lineage}_{segment}_{resolution}_meta.json"

rule files:
    params:
        input_fasta = path_to_fauna+"/{lineage}_{segment}.fasta",
        # dropped_strains = "config/dropped_strains_{lineage}.txt",
        reference = "config/{lineage}_{segment}_outgroup.gb",
        colors = "config/colors.tsv",
        auspice_config = "config/auspice_config.json"

files = rules.files.params


rule parse:
    message: "Parsing fasta into sequences and metadata"
    input:
        sequences = files.input_fasta
    output:
        sequences = "results/sequences_{lineage}_{segment}.fasta",
        metadata = "results/metadata_{lineage}_{segment}.tsv"
    params:
        fasta_fields =  "strain virus isolate_id date region country division passage authors age gender"

    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields}
        """


rule select_strains:
    input:
        metadata = lambda w:expand("results/metadata_{lineage}_{segment}.tsv", segment=segments, lineage=w.lineage)
    output:
        strains = "results/strains_{build}_{lineage}_{resolution}.txt",
    params:
        viruses_per_month = vpm
    shell:
        """
        python scripts/prepare.py --metadata {input.metadata} \
                                  --segments {segments} \
                                  --resolution {wildcards.resolution} --lineage {wildcards.lineage} \
                                  --viruses_per_month {params.viruses_per_month} \
                                  --output {output.strains}
        """


rule filter:
    input:
        metadata = rules.parse.output.metadata,
        sequences = 'results/sequences_{lineage}_{segment}.fasta',
        strains = rules.select_strains.output.strains
    output:
        sequences = 'results/sequences_{build}_{lineage}_{segment}_{resolution}.fasta'
    run:
        from Bio import SeqIO
        with open(input.strains) as infile:
            strains = set(map(lambda x:x.strip(), infile.readlines()))
        with open(output.sequences, 'w') as outfile:
            for seq in SeqIO.parse(input.sequences, 'fasta'):
                if seq.name in strains:
                    SeqIO.write(seq, outfile, 'fasta')



rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned_{build}_{lineage}_{segment}_{resolution}.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/treeraw_{build}_{lineage}_{segment}_{resolution}.nwk"
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
        tree = "results/tree_{build}_{lineage}_{segment}_{resolution}.nwk",
        node_data = "results/branchlengths_{build}_{lineage}_{segment}_{resolution}.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4
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
        node_data = "results/ntmuts_{build}_{lineage}_{segment}_{resolution}.json"
    params:
        inference = "joint"
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
        node_data = "results/aamuts_{build}_{lineage}_{segment}_{resolution}.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """

rule titers:
    input:
        tree = rules.refine.output.tree,
        titers = titer_data
    output:
        tree_model = "results/HITreeModel_{build}_{lineage}_{segment}_{resolution}.json",
        #subs_model = "{subtype}/results/HI_subs_model.json"
    shell:
        """
        augur titers --tree {input.tree}\
            --titers {input.titers}\
            --titer-model tree \
            --output {output.tree_model}
        """

rule export:
    input:
        tree = rules.refine.output.tree,
        node_data = rules.refine.output.node_data,
        metadata = rules.parse.output.metadata,
        nt_muts = rules.ancestral.output,
        aa_muts = rules.translate.output,
        tree_model = rules.titers.output.tree_model,
        auspice_config = files.auspice_config
    output:
        auspice_tree = "auspice/flu_{build}_{lineage}_{segment}_{resolution}_tree.json",
        auspice_meta = "auspice/flu_{build}_{lineage}_{segment}_{resolution}_meta.json"
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} {input.nt_muts} {input.aa_muts} {input.tree_model}\
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_meta}
        """
