
"""
Developed for python3.8
justingreeblatt@github | last updated 26/06/2022

This iterates over the rna records in a single file and calls blast (rev_blast.py) on each of them.

      python3 rnaWalk.py {fasta of rna from genomic input} {outFile path} {control outfile path}

dependacies are :
python libraries:
biopython==1.78

CLI:
ncbiblast+
"""

from subprocess import Popen, PIPE
import os
import sys
from Bio import SeqIO
import re
from settings.directories import REV_BLAST_PATH
from settings.logs import RNA_WALK_LOG_PATH, RNA_WALK_LOG_LEVEL
import logging

logging.basicConfig(filename = RNA_WALK_LOG_PATH, level = RNA_WALK_LOG_LEVEL)

#MRNA
IN_FILE_MRNA = sys.argv[1]

#outFiles
OUT_FILE = sys.argv[2]
CONTROL_OUT_FILE = sys.argv[3]

logging.debug("STARTING rnaWalk: mrnaIn hitsOut controlHitsOut \t{}\t{}\t{}\t{}".format(sys.argv[1],sys.argv[2],sys.argv[3],time()))

#create output files
new = open(OUT_FILE,'x')
newControl = open(CONTROL_OUT_FILE, 'x')
new.write("geneId, geneStart, geneEnd, geneStrand, exons, queryStart, queryEnd, subjectStart, subjectEnd, matchlength, matchPct\n")
newControl.write("geneId, geneStart, geneEnd, exons, geneStrand, queryStart, queryEnd, subjectStart, subjectEnd, matchlength, matchPct\n")
new.close()
newControl.close()

rnaHandler = open(IN_FILE_MRNA)

for rna in SeqIO.parse(rnaHandler, "fasta"):
        
    #Getting data of the gene
    m = re.search(r"^(?P<name>.*?)\s.*join\((?P<location>.*?)\)\)\].*", rna.description)
    geneId = m.group("name")
    exons =  m.group("location").replace(',','-')
    geneStart = m.group("location").split("..")[0]
    geneEnd = m.group("location").split("..")[-1]
    geneStrand = ''
    if geneStart < geneEnd:
        geneStrand = '-'
    else:
        geneStrand = '+'

    geneSeq = str(rna.seq)

    #Writing gene data to temporary file 
    tempFilename = "temp_fasta_" + geneId.split('|')[-1] + '.fa'
    tempGeneFasta = open(tempFilename, 'w')
    fastaHeader = ">" + ",".join([geneID, str(geneStart), str(geneEnd), str(geneStrand), exons]) + "\n"
    tempGeneFasta.write(fastaHeader)
    tempGeneFasta.write(geneSeq)
    tempGeneFasta.close()

    #Running revBlast on gene
        
    command = ["sudo", "python3", REV_BLAST_PATH, tempFilename, OUT_FILE, "minus"]
    p = Popen(command, stdout = PIPE)
    p.wait()
        
    #getControl

    controlCommand = ["sudo", "python3", REV_BLAST_PATH, tempFilename, CONTROL_OUT_FILE, "plus"]
    pc = Popen(controlCommand, stdout = PIPE)
    pc.wait()

    #Removing temporary file
    os.remove(tempFilename)

logging.debug("ENDING rnaWalk: genomeIn gtfIn hitsOut controlHitsOut time\t{}\t{}\t{}\t{}".format(sys.argv[1], sys.argv[2], sys.argv[3] ,sys.argv[4] ,time()))
rnaHandler.close()
