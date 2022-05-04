from os import environ
from os.path import join
from time import time

APLICATION_PATH = join("/home/greenblattcloud2/", "blastWeb")
DATA_PATH = join("/home/greenblattcloud2/", "data")



GTF_FOLDER = join(DATA_PATH, "gtf")
GENOME_FOLDER = join(DATA_PATH, "genomes")
GENOME_WALK_PATH = join(DATA_PATH, "bin", "genomeWalk.py")
GENOME_WALK_FOLDER = join(DATA_PATH,"genomeWalk")
REPEAT_MASK_FOLDER = join(DATA_PATH, "repeatMask")
HISTOGRAMS_FOLDER = join(DATA_PATH, "histograms")
HISTOGRAMS2D_FOLDER = join(DATA_PATH, "histograms2D")


REV_BLAST_PATH = join(APLICATION_PATH, "bin", "revBlast.py")
BLAST_REV_TEMP_DIR = join(DATA_PATH,"temp")

ENSEMBL_HTML_PATH = join(DATA_PATH, "static", "ensemblGenomes.html")
FILE_DIRECTORIES = (DATA_PATH, "static", "localDataDirectories.json")
LOG_FOLDER = join(APLICATION_PATH, "logs")

