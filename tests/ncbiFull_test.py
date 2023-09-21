import sys, os, json
sys.path.insert(1, os.path.join(os.environ.get("HOME"), "blastWeb/bin"))
from ncbiData import ncbiData
import logging

class TestNcbiData:

    platypus = ncbiData("Anas_platyrhynchos-GCF_015476345.1_ZJU1.0.json")

    def test_full(self):

        TestNcbiData.platypus.runMrnaBlast()
        TestNcbiData.platypus.runPremrnaBlast()
        TestNcbiData.platypus.generateHistograms()

if __name__ == "__main__":
    T = TestNcbiData()
    T.test_full()