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

    def test_ncbiDataInit(self):
        assert TestNcbiData.platypus.id == TestNcbiData.assemblieName.replace('/','_')

    def getAndDelete(self, key):

        TestNcbiData.platypus.getResource(key)
        keyFileDir = TestNcbiData.platypus.fileDirectories[key]
        assert os.path.isfile(keyFileDir)
        TestNcbiData.platypus.deleteResource(key)
        assert not os.path.isfile(keyFileDir)

    def test_gtf(self):
        self.getAndDelete("gtf")
    def test_genome(self):
        self.getAndDelete("dna")
    def test_mrna(self):
        self.getAndDelete("mrna")
    def test_GenomeWalk(self):
        TestNcbiData.platypus.runGenomeWalk()
        assert os.path.isfile(TestNcbiData.platypus.fileDirectories["genomeWalk"])
        assert os.path.isfile(TestNcbiData.platypus.fileDirectories["genomeWalkControl"])
    def test_rnaWalk(self):
        TestNcbiData.platypus.runRnaWalk()
        assert os.path.isfile(TestNcbiData.platypus.fileDirectories["rnaWalk"])
        assert os.path.isfile(TestNcbiData.platypus.fileDirectories["rnaWalkControl"])

if __name__ == "__main__":
    t = TestNcbiData()
    g = t.platypus
    g.getResource("dna")
    g.getResource("gtf")
