import sys, os, json
sys.path.insert(1, os.path.join(os.environ.get("HOME"), "blastWeb/bin"))
from ncbiData import ncbiData
from settings.directories import NCBI_DUMMY_DATA
import logging
from settings import dConfig

class TestNcbiData:

    
    dummyData = json.load(open(NCBI_DUMMY_DATA))
    assemblieName = list(dummyData)[0]
    assemblieLinks = dummyData[assemblieName]
    platypus = ncbiData(assemblieName, assemblieLinks)

    def test_full(self):

        #TestNcbiData.platypus.runMrnaBlast()
        TestNcbiData.platypus.runPremrnaBlast()

if __name__ == "__main__":
    T = TestNcbiData()
    T.test_full()

