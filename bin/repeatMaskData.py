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


#Download stuff from the internet into the current working directory
def getAllRepeatMaskLinks(url = REPEAT_MASK_BASE_URL):
    raw = get(url).text
    return findall(GET_REPEAT_MASK_SPECIES_REGEX ,raw)

class repeatMaskData():

    def __init__(self,url):
        raw = get(url).text
        self.genome, self.genomeDate =  f = findall(r'h3>(.*?) - (.*?) - ',raw)[0]
        self.species = findall(r'<TITLE>.*?\[ (.*?) \].*?<',raw)[0]
        self.repeatMaskerOutLink = findall(r'href=\"(.*?\.fa\.out\.gz)\"',raw)[0]

        if not self.repeatMaskerOutLink.startswith("https://www.repeatmasker.org"):
            self.repeatMaskerOutLink = "https://www.repeatmasker.org" + self.repeatMaskerOutLink

        self.initiated = False

    def getRepeatData(self):

        cwd = getcwd()
        chdir(REPEAT_MASK_FOLDER)
        filename = downloadFromURL(self.repeatMaskerOutLink) 
        newFilename = "/".join(filename.split('/')[:-1]) + self.species + ".out.gz"
        newFilename.replace(" ", "_")
        rename(filename, newFilename)
        chdir(cwd)


#        self.familyDist = Counter()
#        self.classDist = Counter()

#        h = gzip.open(join(REPEAT_MASK_FOLDER, filename))
#        self.sizeDist = []

#        self.outDict = {}

#        for n,l in enumerate(h):
#            if n > 2:
#                l = str(l)
#                l = l.strip("b\'").strip("\\n").split()
#                if l[-5] in self.outDict:
#                    self.outDict[l[-5]].append((int(l[6]) - int(l[5])))
#                else:
#                    self.outDict[l[-5]] = []

#        self.sums = {a:mean(self.outDict[a]) for a in self.outDict}
#        self.repeatDist = [a for b in self.outDict for a in self.outDict[b]]

#        h.close()
#        self.repeatData = None

#        self.initated = True



def getAllData():
    links = getAllRepeatMaskLinks()
    species = []
    for l in links:
        organism = repeatMaskData(l)
        print(organism.species)
        species.append(organism)
    return species
