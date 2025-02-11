configfile: "config.yaml"

def get_texts(wildcards):
    return config["texts"][wildcards.text]

# REQUIRED: Run the download rule first!
rule all:
    input:
        expand("Graphs/char_sorted/{text}.svg", text=config["texts"]),
        expand("Graphs/magnitude_sorted/{text}.svg", text=config["texts"]),
        "Graphs/combined_chars/all_char.svg",
        "Graphs/combined_chars/all_magnitude.svg"

rule text_to_chars:
    input:
        get_texts
    output:
        "CharDataOutput/{text}_chars.txt"
    log:
        out = 'Logs/text_to_chars/{text}_stdout.log',
        err = 'Logs/text_to_chars/{text}_stderr.err'
    # group: "individual_texts"      # used in cluster computing, like-grouped rules are submitted to the same node (see https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html for DAG info)
    threads: 4  
    conda:
        "environments/char_frequency_env.yaml"
    shell:
        "bash Scripts/text_to_chars.sh {input} {output} 2> {log.err} 1> {log.out}"

rule combine_char_data_to_file:
    input:
        expand("CharDataOutput/{text}_chars.txt", text=config["texts"])
    output:
        "CombinedCharOutput/all_text_chars.txt"
    log:
        out = 'Logs/combine_char_data_to_file/stdout.log',
        err = 'Logs/combine_char_data_to_file/stderr.err'
    # group: "combined_texts"
    threads: 5
    conda:
        "environments/char_frequency_env.yaml"
    shell:
        "bash Scripts/combine_char_data.sh {output} {input} 2> {log.err} 1> {log.out}"

rule bar_graph_by_char:
    input:
        "CharDataOutput/{text}_chars.txt"
    output:
        "Graphs/char_sorted/{text}.svg"
    params:
        mode="char"
    log:
        out = 'Logs/bar_graph_by_char/{text}_stdout.log',
        err = 'Logs/bar_graph_by_char/{text}_stderr.err'
    # group: "individual_texts"
    threads: 2
    conda:
        "environments/char_frequency_env.yaml"
    script:
        "Scripts/SortAndGraph.py"

rule bar_graph_by_magnitude:
    input:
        "CharDataOutput/{text}_chars.txt"
    output:
        "Graphs/magnitude_sorted/{text}.svg"
    params:
        mode="magnitude"
    log:
        out = 'Logs/bar_graph_by_magnitude/{text}_stdout.log',
        err = 'Logs/bar_graph_by_magnitude/{text}_stderr.err'
    # group: "individual_texts"
    threads: 2
    conda:
        "environments/char_frequency_env.yaml"
    script:
        "Scripts/SortAndGraph.py"

rule combined_graph_char:
    input:
        "CombinedCharOutput/all_text_chars.txt"
    output:
        "Graphs/combined_chars/all_char.svg"
    params:
        mode="char"
    log:
        out = 'Logs/combined_graph_char/stdout.log',
        err = 'Logs/combined_graph_char/stderr.err'
    # group: "combined_texts"
    threads: 2
    conda:
        "environments/char_frequency_env.yaml"
    script:
        "Scripts/SortAndGraph.py"

rule combined_graph_magnitude:
    input:
        "CombinedCharOutput/all_text_chars.txt"
    output:
        "Graphs/combined_chars/all_magnitude.svg"
    params:
        mode="magnitude"
    log:
        out = 'Logs/combined_graph_magnitude/stdout.log',
        err = 'Logs/combined_graph_magnitude/stderr.err'
    # group: "combined_texts"
    threads: 2
    conda:
        "environments/char_frequency_env.yaml"
    script:
        "Scripts/SortAndGraph.py"
