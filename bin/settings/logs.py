from os.path import join
from settings.directories import LOG_FOLDER

#MASTER SWITCH FOR LOGGING LEVEL
MASTER_LOG_LEVEL = "debug"

#GENERATING UNIQUE LOG_ID FROM LIST OF BRAZILIAN CITYS
LOG_ID_TEMPLATE = join(LOG_FOLDER, "logIDs","cityIDs.csv")
LOG_ID_LIST_PATH = join(LOG_FOLDER, "logIDs","cityIDsNotUsed.csv")
loadedIDs = open(LOG_ID_LIST_PATH).read().split(',')
if len(loadedIDs) < 5:
    loadedIDs = open(LOG_ID_TEMPLATE).read().split(',')
LOG_ID = choice(loadedIDs)
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
