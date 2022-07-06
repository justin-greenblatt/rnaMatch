"""
justingreenblatt@github.com |last updated 9/5/2022
Logging for this project has twist. Each run of the code creates a folder insides directories.LOG_FOLDER 
with a name from randomly picked city of Brazil. Path to log and log LEVEL is specified here. Default
setting is set by MASTER_LOG_LEVEL.
"""
from os.path import join
from settings.directories import LOG_FOLDER
from os import mkdir
import logging
from random import choice
from datetime import datetime
from time import strftime

#MASTER SWITCH FOR LOGGING LEVEL
MASTER_LOG_LEVEL = logging.DEBUG

#GENERATING UNIQUE LOG_ID FROM TIMESTAMP
t = datetime.now()
LOG_ID = "blastWeb_{}D{}M{}Y_{}h{}m{}s".format(t.day, t.month, t.year, t.hour, t.minute, t.second)
mkdir(join(LOG_FOLDER, LOG_ID))
L = lambda x: join(LOG_FOLDER, LOG_ID, "{}_{}.txt".format(LOG_ID, x))
TIME = lambda: strftime("%c")

#DEFINING LOGS AND LEVELS
MAIN_LOG_PATH = L("main")
MAIN_LOG_LEVEL = MASTER_LOG_LEVEL

GENOME_WALK_LOG_PATH = L("gw")
GENOME_WALK_LOG_LEVEL = MASTER_LOG_LEVEL 

RNA_WALK_LOG_PATH = L("rw")
RNA_WALK_LOG_LEVEL = MASTER_LOG_LEVEL

HISTOGRAM_LOG_PATH = L("hist")
HISTOGRAM_LOG_LEVEL = MASTER_LOG_LEVEL

HISTOGRAM2D_LOG_PATH = L("hist2d")
HISTOGRAM2D_LOG_LEVEL = MASTER_LOG_LEVEL

ENSEMBL_GENOME_LOG_PATH = L("ensembl")
ENSEMBL_GENOME_LOG_LEVEL = MASTER_LOG_LEVEL

NCBI_GENOME_LOG_PATH = L("ncbi")
NCBI_GENOME_LOG_LEVEL = MASTER_LOG_LEVEL

NCBI_SCRAPE_LOG_PATH = L("ncbi")
NCBI_SCRAPE_LOG_LEVEL = MASTER_LOG_LEVEL

REPEAT_MASK_LOG_PATH = L("repeatMask")
REPEAT_MASK_LOG_LEVEL = MASTER_LOG_LEVEL

REV_BLAST_LOG_PATH = L("repeatMask")
REV_BLAST_LOG_LEVEL = MASTER_LOG_LEVEL

DOWNLOADS_LOG_PATH = L("downloads")
DOWNLOADS_LOG_LEVEL = MASTER_LOG_LEVEL
