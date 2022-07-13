from subprocess import Popen, PIPE
import os
from time import time
import sys
import logging
import logging.config
import Bio.Blast.NCBIXML as BlastReader
from settings.directories import LOGGING_CONFIG, PROCESSES_CONFIG, DIRECTORIES_CONFIG
from configparser import ConfigParser, ExtendedInterpolation

pConfig = ConfigParser()
pConfig.read(PROCESSES_CONFIG)

dConfig = ConfigParser(interpolation = ExtendedInterpolation())
dConfig.read(DIRECTORIES_CONFIG)

logging.config.fileConfig(LOGGING_CONFIG)
logger = logging.getLogger("revBlast")
logger.debug(f'Running revBlast {" ".join(sys.argv)}')

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
    logger.error("Entry error in BlastRev. Exiting")
    sys.exit(0)

name = tempPrefix + sys.argv[1].split("/")[-1].split(".")[0]
if not os.path.exists(IN_FILE):
    sys.exit(0)
blastResultsFile = name + '.' + pConfig["blastn"]["OUT_FILE_SUFFIX"]
blastResultsPath = os.path.join(dConfig["revBlast"]["REV_BLAST_TEMP_FOLDER"], blastResultsFile)
#CREATING BLAST DATABASE
makeDbCommand = list([a for a in 

        [pConfig["makeblastdb"]["USER"],
        pConfig["makeblastdb"]["USER_FLAGS"],
        pConfig["makeblastdb"]["RUN_PATH"],
        "-in" , IN_FILE, 
        "-out", name,
        "-dbtype", pConfig["makeblastdb"]["DB_TYPE"],
        pConfig["makeblastdb"]["EXTRA_FLAGS"],
        pConfig["makeblastdb"]["EXTRA_FLAG_VALUE_PAIRS"]]
        if a])

p = Popen(makeDbCommand, stdout = PIPE)
logger.debug(f'Creating Blast {" ".join(makeDbCommand)}')
p.wait()

#RUNING BLAST ALIGNER
blastnCommand = list([a for a in 

        [pConfig["blastn"]["USER"],
        pConfig["blastn"]["USER_FLAGS"],
        pConfig["blastn"]["RUN_PATH"],
        "-query",IN_FILE, 
        "-strand",STRAND,
        "-db", name,
        "-out", blastResultsPath,
        "-outfmt", pConfig["blastn"]["OUT_FORMAT"],
        pConfig["blastn"]["EXTRA_FLAGS"],
        pConfig["blastn"]["EXTRA_FLAG_VALUE_PAIRS"]]
        if a])

logger.debug(f'Running blastn {" ".join(blastnCommand)}')
b = Popen(blastnCommand, stdout = PIPE)
b.wait()

#WRITING OUTPUT
fastaHandler = open(IN_FILE)
fastaHeader = fastaHandler.readline()
fastaHandler.close()
summary = ""
blast = list(BlastReader.parse(open(blastResultsPath)))
blastRecords = []
try:
    blastRecords = blast[0].alignments[0].hsps
except:
    logger.error(f'Not able to read blast results from {blastResultsPath}')
    sys.exit(0)

logger.debug(f'Writing out blast results: nHits: {len(blastRecords)}')
for a in blastRecords:
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
        logger.debug(f"removing file: {f}")
