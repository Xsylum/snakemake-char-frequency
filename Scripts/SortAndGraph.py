import matplotlib.pyplot as plt
import numpy as np
from snakemake.script import snakemake
import sys

# Functions
def charSort(characters, occurrences):
    records=[] # (char, # of occurrences) pairs
    for record in stringRecords:
        recordValues=record.split()
        records+=[(recordValues[0], recordValues[1])]

    records.sort()
    for record in records:
        characters+=[record[0]]
        occurrences+=[int(record[1])]

def magnitudeSort(characters, occurrences):
    records=[]  # (char, # of occurrences) pairs
    for record in stringRecords:
        recordValues=record.split()
        records+=[(int(recordValues[1]), recordValues[0])]

    records.sort()
    records.reverse()
    for record in records:
        characters+=[record[1]]
        occurrences+=[int(record[0])]

# global variables
characters=[]
occurrences=[]

snakeInput=snakemake.input[0]
sortFormat=snakemake.params.mode  # Either "char" or "magnitude"

# # MAIN
# with open(snakemake.log[0], 'w') as logOut, open(snakemake.log[1], 'w') as logErr:
#     # preparing log files
#     sys.stdout = logOut
#     sys.stderr = logErr

# getting input
with open(snakeInput, 'r') as file:
    stringRecords = file.readlines()

# graph title prep
tmpArray=snakeInput.split("/")
tmpArray=tmpArray[1].split("_")
graphTitle=tmpArray[0] + " Character Contents"

# sort input and complete title
match sortFormat:
    case "char":
        charSort(characters, occurrences)
        graphTitle += " (Sorted by Char)"
    case "magnitude":
        magnitudeSort(characters, occurrences)
        graphTitle += " (Sorted by most occurrences)"

# complete bar graph and save it to output
plt.title(graphTitle)
plt.bar(characters, occurrences)
plt.ylabel("occurrences in text")
plt.xlabel("character")
plt.ylim(0, max(occurrences) + 1)   # bottom of graph at 0, whitespace above max # of occurrences
plt.savefig(snakemake.output[0], dpi=1000)