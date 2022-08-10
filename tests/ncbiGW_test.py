import sys, os, json
sys.path.insert(1, os.path.join(os.environ.get("HOME"), "blastWeb/bin"))
from ncbiData import ncbiData
import pytest
from settings.directories import NCBI_DUMMY_DATA
import logging

class TestNcbiData:

    
    dummyData = json.load(open(NCBI_DUMMY_DATA))
    assemblieName = list(dummyData)[0]
    assemblieLinks = dummyData[assemblieName]
    platypus = ncbiData(assemblieName, assemblieLinks)

    def test_ncbiData(self):
        TestNcbiData.platypus.runGenomeWalk()
       
