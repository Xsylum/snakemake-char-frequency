configfile: "workflow/snakefile_config.yaml"

def get_texts(wildcards):
    return config["texts"][wildcards.text]

# REQUIRED: Run the download rule first!
rule all:
    input:
        expand("workflow/data/graphs/char_sorted/{text}.svg", text=config["texts"]),
        expand("workflow/data/graphs/magnitude_sorted/{text}.svg", text=config["texts"]),
        "workflow/data/graphs/combined_chars/all_char.svg",
        "workflow/data/graphs/combined_chars/all_magnitude.svg"

rule text_to_chars:
    input:
        get_texts
    output:
        "workflow/data/char_data_output/{text}_chars.txt"
    log:
        out = 'workflow/logs/snakemake_logs/text_to_chars/{text}_stdout.log',
        err = 'workflow/logs/snakemake_logs/text_to_chars/{text}_stderr.err'
    # group: "individual_texts"      # used in cluster computing, like-grouped rules are submitted to the same node (see https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html for DAG info)
    threads: 4  
    conda:
        "environments/char_frequency_env.yaml"
    shell:
        "bash workflow/scripts/text_to_chars.sh {input} {output} 2> {log.err} 1> {log.out}"

rule combine_char_data_to_file:
    input:
        expand("workflow/data/char_data_output/{text}_chars.txt", text=config["texts"])
    output:
        "workflow/data/combined_char_output/all_text_chars.txt"
    log:
        out = 'workflow/logs/snakemake_logs/combine_char_data_to_file/stdout.log',
        err = 'workflow/logs/snakemake_logs/combine_char_data_to_file/stderr.err'
    # group: "combined_texts"
    threads: 5
    conda:
        "environments/char_frequency_env.yaml"
    shell:
        "bash workflow/scripts/combine_char_data.sh {output} {input} 2> {log.err} 1> {log.out}"

rule bar_graph_by_char:
    input:
        "workflow/data/char_data_output/{text}_chars.txt"
    output:
        "workflow/data/graphs/char_sorted/{text}.svg"
    params:
        mode="char"
    log:
        out = 'workflow/logs/snakemake_logs/bar_graph_by_char/{text}_stdout.log',
        err = 'workflow/logs/snakemake_logs/bar_graph_by_char/{text}_stderr.err'
    # group: "individual_texts"
    threads: 2
    conda:
        "environments/char_frequency_env.yaml"
    script:
        "workflow/scripts/SortAndGraph.py"

rule bar_graph_by_magnitude:
    input:
        "workflow/data/char_data_output/{text}_chars.txt"
    output:
        "workflow/data/graphs/magnitude_sorted/{text}.svg"
    params:
        mode="magnitude"
    log:
        out = 'workflow/logs/snakemake_logs/bar_graph_by_magnitude/{text}_stdout.log',
        err = 'workflow/logs/snakemake_logs/bar_graph_by_magnitude/{text}_stderr.err'
    # group: "individual_texts"
    threads: 2
    conda:
        "environments/char_frequency_env.yaml"
    script:
        "workflow/scripts/SortAndGraph.py"

rule combined_graph_char:
    input:
        "workflow/data/combined_char_output/all_text_chars.txt"
    output:
        "workflow/data/graphs/combined_chars/all_char.svg"
    params:
        mode="char"
    log:
        out = 'workflow/logs/snakemake_logs/combined_graph_char/stdout.log',
        err = 'workflow/logs/snakemake_logs/combined_graph_char/stderr.err'
    # group: "combined_texts"
    threads: 2
    conda:
        "environments/char_frequency_env.yaml"
    script:
        "workflow/scripts/SortAndGraph.py"

rule combined_graph_magnitude:
    input:
        "workflow/data/combined_char_output/all_text_chars.txt"
    output:
        "workflow/data/graphs/combined_chars/all_magnitude.svg"
    params:
        mode="magnitude"
    log:
        out = 'workflow/logs/snakemake_logs/combined_graph_magnitude/stdout.log',
        err = 'workflow/logs/snakemake_logs/combined_graph_magnitude/stderr.err'
    # group: "combined_texts"
    threads: 2
    conda:
        "environments/char_frequency_env.yaml"
    script:
        "workflow/scripts/SortAndGraph.py"
