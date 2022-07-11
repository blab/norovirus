rule all:
    input:
        auspice_json = "auspice/norovirus.json",

reference = "config/norovirus_outgroup.gb",
auspice_config = "config/auspice_config.json"

rule parse:
    input:
        sequences = "data/sequences_vipr.fasta",
    output:
        sequences = "results/sequences.fasta",
        metadata = "results/metadata_vipr.tsv"
    params:
        fields = ["strain","strain_name","segment","date","host","country"]
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fix-dates monthfirst \
            --fields {params.fields}
        """
rule clean_dates:
    input:
        metadata = "results/metadata_vipr.tsv"
    output:
        result = "results/metadata_vipr_cleaned_dates.tsv"
    params:
        columns = ["strain","strain_name","segment","date","host","country"]
    shell:
        """
        ./scripts/tsv-to-ndjson < {input.metadata} \
            | ./scripts/transform-date-fields --expected-date-formats "%Y_%m_%d" "%Y_%m_%dT%H:%M:%SZ" "%Y_%m" "%Y-%m-%d" "%Y-XX-XX" "%Y-%m" --date-fields date \
            | ./scripts/ndjson-to-tsv --metadata-columns {params.columns} --metadata {output.result}
        """
rule prepare_genotype_metadata:
    input:
        metadata = expand("data/genomicdetective_results{i}.csv", i = [1,2,3])
    output:
        result = "results/metadata_genomicdetective.tsv"
    params:
        fields = "ORF2_type,strain"
    shell:
        """
        csvtk concat {input.metadata} \
            | csvtk rename -f "ORF2 type" -n "ORF2_type" \
            | csvtk cut -f {params.fields} \
            | csvtk --out-tabs replace -f "strain" -p "\|.*$" -r "" > {output.result}
        """

rule join_metadata:
    input:
        genomicdetective_metadata = "results/metadata_genomicdetective.tsv",
        vipr_metadata = "results/metadata_vipr_cleaned_dates.tsv"
    output:
        result = "results/metadata.tsv"
    shell:
        """
        csvtk -t join -f "strain" {input.vipr_metadata} {input.genomicdetective_metadata} --left-join --na "NA" > {output.result}
        """

rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - from {params.min_date} onwards
          - excluding strains in {input.exclude}
          - minimum genome length of {params.min_length} (67% of Norovirus virus genome)
        """
    input:
        sequences = "results/sequences.fasta",
        metadata = "results/metadata.tsv",
        exclude = "config/dropped_strains.txt"
    output:
        sequences = "results/filtered.fasta",
        metadata = "results/filtered_metadata.tsv"
    params:
        group_by = "year ORF2_type",
        sequences_per_group = 30,
        min_date = 1950,
        min_length = 5032
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --query "ORF2_type in ['GII.6', 'GII.4', 'GII.2', 'GII.3', 'GII.17']" \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --min-date {params.min_date} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-length {params.min_length}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = reference
    output:
        alignment = "results/aligned.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --nthreads 4
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads 4
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
        metadata = "results/filtered_metadata.tsv"
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --root "best" \
            --coalescent {params.coalescent} \
            --clock-filter-iqd {params.clock_filter_iqd} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --date-confidence \
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt_muts.json"
    params:
        inference = "joint"
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
        reference = reference
    output:
        node_data = "results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
        """

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = "results/filtered_metadata.tsv",
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        auspice_config = auspice_config
    output:
        auspice_json = rules.all.input.auspice_json,
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """
