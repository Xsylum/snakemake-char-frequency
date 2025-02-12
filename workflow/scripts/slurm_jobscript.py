#!/usr/bin/env python3
"""
MIT License

Copyright (c) His Majesty The King in Right of Canada, as represented by the Minister of Agriculture and Agri-Food, 2024 | Copyright (c) Sa Majesté le Roi du chef du Canada, représentée par le ministre de l’Agriculture et de l’Agroalimentaire, 2024

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


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

logging.basicConfig(format="%(asctime)s - %(message)s",
                    datefmt="%d-%b-%y %H:%M:%S",
                    level=logging.INFO,
                    filename="workflow/logs/slurm_logs/slurm_jobscript.log",
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

# haven't figured out how to get this work in the config.yaml file
# if "cluster_cancel" in job_properties["resources"]:
#     cluster_cancel = f"--parsable --cluster-cancel {job_properties["resources"]["cluster_cancel"]}"
# else:
#     cluster_cancel = ""
cluster_cancel = ""

cmdline = f"""sbatch -N {cluster_nodes} -n {cluster_ntasks} -p {cluster_partition} \
                     -A {cluster_account} -c {threads} -t {resources_runtime} \
                     -J {rulename}.{jobid} \
                     -o workflow/logs/slurm_logs/{rulename}.{jobid}.{timestamp}.log \
                     -e workflow/logs/slurm_logs/{rulename}.{jobid}.{timestamp}.err \
                     {cluster_cancel} {resources_mem} {clusters} {qos} {jobscript}"""

# add the following to 2nd-to-last line of cmdline for email notifications: --mail-user={cluster_mail_user} --mail-type={cluster_mail_type} \

print(cmdline)
logger.info(f"rule: {rulename}")
logger.info(f"job_properties['resources']: {job_properties['resources']}")
logger.info(f"command line: {cmdline}")
logger.info("test")
os.system(cmdline)
