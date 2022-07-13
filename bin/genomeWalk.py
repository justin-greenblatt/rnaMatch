
"""
Developed for python3.9
justingreeblatt@github | last updated 12/06/2022

This iterates over the genome and annotation and calls blast (rev_blast.py) on each gtf coordinate.

Its function is mainly to iterate in a paired manner over the genome and coordinates. IT IS ASSUMED THAT THE
GTF COORDINATES ARE SORTED ACORDING TO THE ORDER OF THE CHROMOSSOMES IN THE GENOME FILE.
This code also manages output of the revBlast calls and groups them in 2 output files for the whole genome.
arg[3] or {OUT_FILE} are the results obtained by using blast of the region on its reverse complement.
arg[4] or {OUT_FILE_CONTROL} has the same content but by using blast of a region on itself.
If you are running this out of the blastWeb repository modify the code or create your own confign settings.
It is a command line tool  / python script that should be used in the following manner

      python3 genomeWalk.py {genome input} {gtf annotation input} {outFile path} {control outfile path}

dependacies are :
python libraries:
biopython==1.78

CLI:
ncbiblast+
"""

from subprocess import Popen, PIPE
from configparser import ConfigParser, ExtendedInterpolation

import os
import sys

from Bio import SeqIO
import re

from settings.directories import LOGGING_CONFIG, PROCESSES_CONFIG, DIRECTORIES_CONFIG

import gzip
import logging
import logging.config

#Setting up Logging
logging.config.fileConfig(LOGGING_CONFIG)
logger = logging.getLogger(__name__)

pConfig = ConfigParser()
pConfig.read(PROCESSES_CONFIG)

dConfig = ConfigParser(interpolation = ExtendedInterpolation())
dConfig.read(DIRECTORIES_CONFIG)

#Genome
IN_FILE_GENOME = sys.argv[1]

#GTF
IN_FILE_GTF = sys.argv[2]

#outFiles
OUT_FILE = sys.argv[3]
CONTROL_OUT_FILE = sys.argv[4]
logger.debug(f"STARTING genomeWalk:{sys.argv[0]} {sys.argv[1]} {sys.argv[2]} {sys.argv[3]} {sys.argv[4]}")

#create output files
new = open(OUT_FILE,'x')
newControl = open(CONTROL_OUT_FILE, 'x')
new.write("geneId,geneStart,geneEnd,geneStrand,queryStart,queryEnd,subjectStart,subjectEnd,matchlength,matchPct\n")
newControl.write("geneId,geneStart,geneEnd,geneStrand,queryStart,queryEnd,subjectStart,subjectEnd,matchlength,matchPct\n")
new.close()
newControl.close()

#open input files as iterators
genomeHandler = None
if IN_FILE_GENOME.endswith("gz"):
    genomeHandler = gzip.open(IN_FILE_GENOME,'rt')
else:
    genomeHandler = open(IN_FILE_GENOME)

genomeIterator = SeqIO.parse(genomeHandler, "fasta")
gtfIterator = None
if IN_FILE_GTF.endswith("gz"):
    gtfIterator = gzip.open(IN_FILE_GTF, 'rt')
else:
    gtfIterator = open(IN_FILE_GTF)

flag = True

chromosome = next(genomeIterator, None)

#function for getting next gene in GTF files.
"""
IF YOU WANT TO FIND OTHER ENTRIES THAT ARE NOT UNDER "gene"
MODIFY THIS CODE.
CHANGE LATER: I would probably make a folder with functions for different annotations and filters
"""

def getNextGene(iterator):
    line = next(iterator)
    while line.startswith('#') or line.split('\t')[2] != "gene":
        line = next(iterator)
        if line == "":
            sys.exit(0)
            break
    return line

gene = getNextGene(gtfIterator)

#Start iterating over genome and annotations
while flag:

    if isinstance(chromosome, type(None)):

        sys.exit(0)
        break

    #skipComments
    while gene.split('\t')[0] == chromosome.id:
        
        #Getting data of the gene
        geneID = re.search(r'gene_id \"(.+?)\"', gene).group(1)
        geneStart = int(gene.split('\t')[3])
        geneEnd = int(gene.split('\t')[4])
        geneStrand = gene.split('\t')[6]
        geneSeq = str(chromosome.seq[max(0, geneStart - int(pConfig["genomeWalk"]["DOWNSTREAM_GENE_FLANK"])):min(geneEnd + int(pConfig["genomeWalk"]["UPSTREAM_GENE_FLANK"]), len(chromosome.seq))])

        #Writing gene data to temporary file 
        tempFilename = "temp_fasta_" + geneID + '.fa'
        tempGeneFasta = open(tempFilename, 'w')
        fastaHeader = ">" + ",".join([geneID, str(geneStart), str(geneEnd), str(geneStrand)]) + "\n"
        tempGeneFasta.write(fastaHeader)
        tempGeneFasta.write(geneSeq)
        tempGeneFasta.close()

        #Running revBlast on gene
        command  = list([a for a in [
                pConfig["revBlast"]["USER"],
                pConfig["revBlast"]["USER_FLAGS"],
                pConfig["revBlast"]["INTERPRETER"],
                pConfig["revBlast"]["INTERPRETER_FLAGS"]]
                        if a])

        testCommand = command + [dConfig["scripts"]["REV_BLAST_PATH"], tempFilename, OUT_FILE, "minus"]
        logger.debug(f"genomeWalk Test running command: {testCommand}")
        p = Popen(testCommand, stdout = PIPE)
        p.wait()
        
        #getControl

        controlCommand = command + [dConfig["scripts"]["REV_BLAST_PATH"], tempFilename, CONTROL_OUT_FILE, "plus"]
        logger.debug(f"genomeWalk Test running command: {controlCommand}")
        pc = Popen(controlCommand, stdout = PIPE)
        pc.wait()

        #Removing temporary file
        os.remove(tempFilename)

        gene = getNextGene(gtfIterator)
    chromosome = next(genomeIterator, None)
logger.debug(f"ENDING genomeWalk:{sys.argv[0]} {sys.argv[1]} {sys.argv[2]} {sys.argv[3]} {sys.argv[4]}")
gtfHandler.close()
genomeHandler.close()

