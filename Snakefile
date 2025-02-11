# https://gcxgce.sharepoint.com/teams/1000645/SitePages/Snakemake-GPSC-Documentation.aspx Under "Add Download Rule/Script"
localrules: download   # run rule on head node, which has access to the internet
include: "MainWorkflow.smk"

# conda:
#     "char_frequency_env.yaml"

# will run on the head node, access to internet
rule download:
    output:
        "conda_packages_downloaded.txt"
    threads: 5
    shell:
        "snakemake --sdm conda --conda-create-envs-only all && touch {output}"