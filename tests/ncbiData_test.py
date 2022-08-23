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
        assert TestNcbiData.platypus.id == TestNcbiData.assemblieName.replace('/','_')
        TestNcbiData.platypus.getResource("dna")
        genomeDir = TestNcbiData.platypus.fileDirectories["dna"]
        assert os.path.isfile(genomeDir)

        TestNcbiData.platypus.getResource("gtf")
        gtfDir = TestNcbiData.platypus.fileDirectories["gtf"]
        assert os.path.isfile(gtfDir)

        TestNcbiData.platypus.getResource("mrna")
        mrnaDir = TestNcbiData.platypus.fileDirectories["mrna"]
        assert os.path.isfile(mrnaDir)

        TestNcbiData.platypus.getResource("gff")
        gffDir = TestNcbiData.platypus.fileDirectories["gff"]
        assert os.path.isfile(gffDir)


if __name__ == "__main__":
    T = TestNcbiData()
    T.test_ncbiData()

