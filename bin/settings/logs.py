"""
justingreenblatt@github.com |last updated 9/5/2022
Logging for this project has twist. Each run of the code creates a folder insides directories.LOG_FOLDER 
with a name from randomly picked city of Brazil. Path to log and log LEVEL is specified here. Default
setting is set by MASTER_LOG_LEVEL.
"""
from os.path import join
from settings.directories import LOG_FOLDER
from os import mkdir

#MASTER SWITCH FOR LOGGING LEVEL
MASTER_LOG_LEVEL = "info"

#GENERATING UNIQUE LOG_ID FROM LIST OF BRAZILIAN CITYS
LOG_ID_TEMPLATE = join(LOG_FOLDER, "logIDs","cityIDs.csv")
LOG_ID_LIST_PATH = join(LOG_FOLDER, "logIDs","cityIDsNotUsed.csv")
loadedIDs = open(LOG_ID_LIST_PATH).read().split(',')
if len(loadedIDs) < 5:
    loadedIDs = open(LOG_ID_TEMPLATE).read().split(',')
LOG_ID = choice(loadedIDs)
mkdir(join(LOG_FOLDER, LOG_ID))
L = lambda x: join(LOG_FOLDER,"{}_{}.txt".format(LOG_ID, x))


#DEFINING LOGS AND LEVELS
MAIN_LOG = L("main")
MAIN_LOG_LEVEL = MASTER_LOG_LEVEL

GENOME_WALK_LOG_PATH = L("gw")
GENOME_WALK_LOG_LEVEL = MASTER_LOG_LEVEL 

HISTOGRAM_LOG_PATH = L("hist")
HISTOGRAM_LOG_LEVEL = MASTER_LOG_LEVEL

HISTOGRAM2D_LOG_PATH = L("hist2d")
HISTOGRAM2D_LOG_LEVEL = MASTER_LOG_LEVEL

ENSEMBL_GENOME_LOG_PATH = L("ensembl")
ENSEMBL_GENOME_LOG_LEVEL = MASTER_LOG_LEVEL

REPEAT_MASK_LOG_PATH = L("repeatMask")
REPEAT_MASK_LOG_LEVEL = MASTER_LOG_LEVEL

REV_BLAST_LOG_PATH = L("repeatMask")
REV_BLAST_LOG_LEVEL = MASTER_LOG_LEVEL
