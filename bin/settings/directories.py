from os import environ
from os.path import join
from time import time

APLICATION_PATH = join(environ.get("HOME"), "rnaMatch")
DATA_PATH = join(environ.get("HOME"), "data")


GTF_FOLDER = join(DATA_PATH, "gtf")
GENOME_FOLDER = join(DATA_PATH, "genomes")
MRNA_FOLDER = join(DATA_PATH, "mrna")
GENOME_WALK_FOLDER = join(DATA_PATH,"genomeWalk")
GENOME_WALK_CONTROL_FOLDER = join(DATA_PATH, "genomeWalkControl")
RNA_WALK_FOLDER = join(DATA_PATH, "rnaWalk")
RNA_WALK_CONTROL_FOLDER = join(DATA_PATH, "rnaWalkControl")
REPEAT_MASK_FOLDER = join(DATA_PATH, "repeatMask")
GFF_FOLDER = join(DATA_PATH, "gff")

HISTOGRAMS_FOLDER = join(DATA_PATH, "histograms")
HISTOGRAMS2D_FOLDER = join(DATA_PATH, "histograms2D")

GENOME_WALK_SIZE_HISTOGRAM_FOLDER = join(HISTOGRAMS_FOLDER, "genomeWalkSize")
GENOME_WALK_CONTROL_SIZE_HISTOGRAM_FOLDER = join(HISTOGRAMS_FOLDER, "genomeWalkSizeControl")
GENOME_WALK_SIZE_HISTOGRAM_2D_FOLDER = join(HISTOGRAMS2D_FOLDER, "genomeWalkSizePct")
GENOME_WALK_CONTROL_SIZE_HISTOGRAM_2D_FOLDER = join(HISTOGRAMS2D_FOLDER, "genomeWalkControlSizePct")

RNA_WALK_SIZE_HISTOGRAM_FOLDER = join(HISTOGRAMS_FOLDER, "rnaWalkSize")
RNA_WALK_CONTROL_SIZE_HISTOGRAM_FOLDER = join(HISTOGRAMS_FOLDER, "rnaWalkControlSize")
RNA_WALK_SIZE_HISTOGRAM_2D_FOLDER = join(HISTOGRAMS2D_FOLDER, "rnaWalkSizePct")
RNA_WALK_SIZE_CONTROL_HISTOGRAM_2D_FOLDER = join(HISTOGRAMS2D_FOLDER, "rnaWalkControlSizePct")

REPEAT_MASK_SIZE_HIST_FOLDER = join(HISTOGRAMS_FOLDER, "repeatMask")

RNA_WALK_PATH = join(DATA_PATH, "bin", "rnaWalk.py")
GENOME_WALK_PATH = join(DATA_PATH, "bin", "genomeWalk.py")
IMAGE_FOLDER = join(APLICATION_PATH, "bin", "static") 
REV_BLAST_PATH = join(APLICATION_PATH, "bin", "revBlast.py")
BLAST_REV_TEMP_DIR = join(DATA_PATH,"temp")
NCBI_RESOURCE_LINKS_PATH = join(APLICATION_PATH, "bin", "ncbiLinks.py")

TEST_FOLDER = join(APLICATION_PATH, "tests")
NCBI_DUMMY_DATA = join(TEST_FOLDER, "dummyNcbiData.json")
TEST_REV_BLAST_IN_FILE = join(TEST_FOLDER, "testGeneData.fa")
TEST_REV_BLAST_OUT_FILE = join(TEST_FOLDER, "testOutMinus.csv")
TEST_REV_BLAST_CONTROL_OUT_FILE = join(TEST_FOLDER, "testOutPlus.csv")

ENSEMBL_HTML_PATH = join(DATA_PATH, "static", "ensemblGenomes.html")
FILE_DIRECTORIES = join(DATA_PATH, "static", "localDataDirectories.json")
LOG_FOLDER = join(APLICATION_PATH, "logs")
ANAGE_DATA_FILE = join(DATA_PATH, "anAge", "anage.txt")
LOGGING_CONFIG = join(APLICATION_PATH, "bin", "settings", "logging.ini")
PROCESSES_CONFIG = join(APLICATION_PATH, "bin", "settings", "processes.ini")
DIRECTORIES_CONFIG = join(APLICATION_PATH, "bin", "settings", "directories.ini")
RESOURCE_FOLDERS = {
            "dna" : GENOME_FOLDER,
            "gtf" : GTF_FOLDER,
            "gff" : GFF_FOLDER,
            "mrna" : MRNA_FOLDER,
            "genomeWalk" : GENOME_WALK_FOLDER,
            "genomeWalkControl" : GENOME_WALK_CONTROL_FOLDER,
            "rnaWalk" : RNA_WALK_FOLDER,
            "rnaWalkControl" : RNA_WALK_CONTROL_FOLDER,
            "repeatMask" : REPEAT_MASK_FOLDER,
            "gwSizeHist" : GENOME_WALK_SIZE_HISTOGRAM_FOLDER,
            "gwCSizeHist" : GENOME_WALK_CONTROL_SIZE_HISTOGRAM_FOLDER,
            "gwSizeHist2d" : GENOME_WALK_SIZE_HISTOGRAM_2D_FOLDER,
            "gwCSizeHist2d" : GENOME_WALK_CONTROL_SIZE_HISTOGRAM_2D_FOLDER,
            "rwSizeHist" : RNA_WALK_SIZE_HISTOGRAM_FOLDER,
            "rwCSizeHist" : RNA_WALK_CONTROL_SIZE_HISTOGRAM_FOLDER,
            "rwSizeHist2d" : RNA_WALK_SIZE_HISTOGRAM_2D_FOLDER,
            "rwCSizeHist2d" : RNA_WALK_SIZE_CONTROL_HISTOGRAM_2D_FOLDER,
            "rmSizeHist" : REPEAT_MASK_SIZE_HIST_FOLDER
            }
