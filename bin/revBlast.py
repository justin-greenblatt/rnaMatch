from subprocess import Popen, PIPE
import os
from time import time
import sys
import logging
import Bio.Blast.NCBIXML as BlastReader
from settings.directories import BLAST_REV_TEMP_DIR
from settings.logs import REPEAT_MASK_LOG_PATH,  REPEAT_MASK_LOG_LEVEL

logging.basicConfig(filename=REV_BLAST_LOG_PATH, level = REV_BLAST_LOG_LEVEL)
logging.debug("Running revBlast: command time [{},{}]".format(" ".join(sys.argv), time()))
#Fasta File
IN_FILE = sys.argv[1]

#Out fasta file
OUT_FILE = sys.argv[2]

#Strand - (plus/minus)
STRAND = sys.argv[3]

tempPrefix = ""
if STRAND == "minus":
    tempPrefix = "blast_temp_"
elif STRAND == "plus":
    tempPrefix = "blast_temp_control_"
else:
    print("Entry error in BlastRev. Exiting")
    sys.exit(0)

name = tempPrefix + sys.argv[1].split("/")[-1].split(".")[0]
if not os.path.exists(IN_FILE):
    sys.exit(0)

#CREATING BLAST DATABASE
print(IN_FILE)
command = ["makeblastdb", "-in" , IN_FILE, "-out", name, "-dbtype", "nucl"]
p = Popen(command, stdout = PIPE)
logging.debug("Creating Blast Database: command time [{}, {}]".format(" ".join(command), time()))
p.wait()

#RUNING BLAST ALIGNER

runBlast = ["blastn", "-query", IN_FILE, "-strand", STRAND, "-db", name,
                "-out", name + ".xml", "-outfmt", "5"]
logging.debug("running blast: command time [{}, {}]".format(" ".join(runBlast), time()))
b = Popen(runBlast, stdout = PIPE)
b.wait()

#WRITING OUTPUT
fastaHandler = open(IN_FILE)
fastaHeader = fastaHandler.readline()
fastaHandler.close()
summary = ""
blast = list(BlastReader.parse(open(name + ".xml")))
blast_records = []
try:
    blast_records = blast[0].alignments[0].hsps
except:
    logging.warning("Not able to read blast results from {}".format(name + ".xml"))
    sys.exit(0)

logging.debug("Writing out blast results: command nHits time [{}, {}, {}]".format(" ".join(runBlast), len(blastRecords) , time()))
for a in blast_records:
    summary += "{},{},{},{},{},{},{}\n".format(fastaHeader.rstrip('\n').lstrip('>'),
                                           a.query_start,
                                           a.query_end,
                                           a.sbjct_start,
                                           a.sbjct_end,
                                           len(a.match),
                                           a.match.count("|")/len(a.match))

summaryOut = open(OUT_FILE, 'a')
summaryOut.write(summary)
summaryOut.close()

for f in os.listdir(os.getcwd()):
    if f.startswith(tempPrefix):
        os.remove(f)
        logging.debug("removing file: relativeFilePath [{}]".format(f))
