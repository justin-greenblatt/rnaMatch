from os import environ
from os.path import join
from time import time

APLICATION_PATH = join("/home/greenblattcloud2/", "blastWeb")
DATA_PATH = join("/home/greenblattcloud2/", "data")


GTF_FOLDER = join(DATA_PATH, "gtf")
GENOME_FOLDER = join(DATA_PATH, "genomes")
MRNA_FOLDER = join(DATA_PATH, "mrna")
GENOME_WALK_FOLDER = join(DATA_PATH,"genomeWalk")
GENOME_WALK_CONTROL_FOLDER = join(DATA_PATH, "genomeWalkControl")
RNA_WALK_FOLDER = join(DATA_PATH, "rnaWalk")
RNA_WALK_CONTROL_FOLDER = join(DATA_PATH, "rnaWalkControl")
REPEAT_MASK_FOLDER = join(DATA_PATH, "repeatMask")

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

ENSEMBL_HTML_PATH = join(DATA_PATH, "static", "ensemblGenomes.html")
FILE_DIRECTORIES = (DATA_PATH, "static", "localDataDirectories.json")
LOG_FOLDER = join(APLICATION_PATH, "logs")
ANAGE_DATA_FILE = join(DATA_PATH, "anAge", "anage.txt")

RESOURCE_FOLDERS = {
            "dna" : GENOME_FOLDER,
            "gtfAnnotation" : GTF_FOLDER,
            "mrna" : MRNA_FOLDER,
            "genomeWalk" : GENOME_WALK_FOLDER,
            "genomeWalkControl" : GENOME_WALK_CONTROL_FOLDER,
            "rnaWalk" : RNA_WALK_FOLDER,
            "rnaControlWalk" : RNA_WALK_CONTROL_FOLDER,
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
