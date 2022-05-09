import gzip
from re import findall
from requests import get
from os import getcwd, chdir, remove, rename
from os.path import join
from collections import Counter
from numpy import mean
from settings.links import REPEAT_MASK_BASE_URL
from settings.directories import REPEAT_MASK_FOLDER
from settings.regularExpressions import GET_REPEAT_MASK_SPECIES_REGEX
from myUtils import downloadFromURL

logging.basicConfig(filename=REPEAT_MASK_LOG_PATH, level = REPEAT_MASK_LOG_LEVEL)

def getAllRepeatMaskLinks(url = REPEAT_MASK_BASE_URL):
    logging.info("Calling getRepeatMaskLinks and fetching all species specific urls from https://www.repeatmasker.org/genomicDatasets/RMGenomicDatasetsAlt.html")
    raw = get(url).text
    return list([repeatMaskData(l) for l in  findall(GET_REPEAT_MASK_SPECIES_REGEX ,raw)])


class repeatMaskData():

    def __init__(self,url):
        logging.info("Initiating repeatMask object from url {}".format(url))
        raw = get(url).text
        self.genome, self.genomeDate =  f = findall(r'h3>(.*?) - (.*?) - ',raw)[0]
        self.species = findall(r'<TITLE>.*?\[ (.*?) \].*?<',raw)[0]
        self.repeatMaskerOutLink = findall(r'href=\"(.*?\.fa\.out\.gz)\"',raw)[0]

        if not self.repeatMaskerOutLink.startswith("https://www.repeatmasker.org"):
            self.repeatMaskerOutLink = "https://www.repeatmasker.org" + self.repeatMaskerOutLink

        self.initiated = False

    def getRepeatData(self):
        logging.info("Fetching repeatMask data for {}|{}".format(self.species, self.genome))
        cwd = getcwd()
        chdir(REPEAT_MASK_FOLDER)
        filename = downloadFromURL(self.repeatMaskerOutLink)
        newFilename = "/".join(filename.split('/')[:-1]) + self.species + ".out.gz"
        newFilename.replace(" ", "_")
        rename(filename, newFilename)
        chdir(cwd)
        self.fileDirectorie = newFilename


def getAllRepeatMaskLinks(url = REPEAT_MASK_BASE_URL):
    logging.info("Calling getRepeatMaskLinks and fetching all species specific urls from https://www.repeatmasker.org/genomicDatasets/RMGenomicDatasetsAlt.html")
    raw = get(url).text
    return list([repeatMaskData(l) for l in  findall(GET_REPEAT_MASK_SPECIES_REGEX ,raw)])
