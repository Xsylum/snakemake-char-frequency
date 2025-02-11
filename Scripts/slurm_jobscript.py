#!/usr/bin/env python3
"""
Author(s): Martin Gauthier
of the Department of Agriculture and Agri-Food, Government of Canada
© His Majesty the King in Right of Canada, as represented by the Minister of Agriculture and Agri-Food Canada 2023
"""


import os
import sys
import datetime
import logging

from snakemake.utils import read_job_properties

logfile_path = os.path.expanduser("~/slurm_jobscript.log")
logging.basicConfig(format="%(asctime)s - %(message)s",
                    datefmt="%d-%b-%y %H:%M:%S",
                    level=logging.INFO,
                    filename=logfile_path,
                    filemode="a")
logger = logging.getLogger('slurm_logger')
logger.setLevel(logging.INFO)
print("testing\n")

jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)

timestamp = datetime.datetime.now().strftime("%Y%m%d.%H%M%S")

threads = job_properties["threads"]
jobid = job_properties["jobid"]
rulename = job_properties["rule"]   # This will be a problem if we desire to submit job groups!
cluster_nodes = job_properties["resources"]["nodes"]
cluster_ntasks = job_properties["resources"]["ntasks"]
cluster_account = job_properties["resources"]["slurm_account"]
cluster_partition = job_properties["resources"]["slurm_partition"]
cluster_mail_type = job_properties["resources"]["mail_type"]
cluster_mail_user = job_properties["resources"]["mail_user"]
resources_runtime = job_properties["resources"]["runtime"]

if "mem_mb_per_cpu" in job_properties["resources"]:
    resources_mem = f'--mem-per-cpu {job_properties["resources"]["mem_mb_per_cpu"]}'
elif "mem_mb" in job_properties["resources"]:
    resources_mem = f'--mem {job_properties["resources"]["mem_mb"]}'
else:
    resources_mem = ""

if "clusters" in job_properties["resources"]:
    clusters = f"--clusters {job_properties['resources']['clusters']}"
else:
    clusters = ""

if "qos" in job_properties["resources"]:
    qos = f"--qos={job_properties['resources']['qos']}"
else:
    qos = ""

cmdline = f"""sbatch -N {cluster_nodes} -n {cluster_ntasks} -p {cluster_partition} \
                     -A {cluster_account} -c {threads} -t {resources_runtime} \
                     -J {rulename}.{jobid} \
                     -o logs/{rulename}.{jobid}.{timestamp}.log \
                     -e logs/{rulename}.{jobid}.{timestamp}.err \
                     --mail-user={cluster_mail_user} --mail-type={cluster_mail_type} \
                     {resources_mem} {clusters} {qos} {jobscript}"""

print(cmdline)
logger.info(f"rule: {rulename}")
logger.info(f"job_properties['resources']: {job_properties['resources']}")
logger.info(f"command line: {cmdline}")
logger.info("test")
os.system(cmdline)
