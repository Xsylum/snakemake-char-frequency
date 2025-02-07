import matplotlib.pyplot as plt
import numpy as np
from snakemake.script import snakemake

characters=[]
occurences=[]

records=[]
snakeInput=snakemake.input[0]
# retrieve character/occurences from file
with open(snakeInput, 'r') as file:
    stringRecords = file.readlines()

tmpArray=snakeInput.split("/")
tmpArray=tmpArray[1].split("_")
graphTitle=tmpArray[0] + " Character Contents (Sorted by most occurences)"
plt.title(graphTitle)

for record in stringRecords:
    recordValues=record.split()
    records+=[(int(recordValues[1]), recordValues[0])]

records.sort()
records.reverse()
for record in records:
    characters+=[record[1]]
    occurences+=[int(record[0])]

plt.bar(characters, occurences)
plt.ylabel("# of occurences in text")
plt.xlabel("character")
plt.ylim(0, max(occurences) + 1)
plt.savefig(snakemake.output[0], dpi=1000)