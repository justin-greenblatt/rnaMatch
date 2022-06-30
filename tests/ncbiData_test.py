import sys, os, json
sys.path.insert(1, os.path.join(os.environ.get("HOME"), "blastWeb/bin"))
from ncbiData import ncbiData
import pytest

def test_ncbiDataInit():
    dummyData = json.load(open("dummyNcbiData.json"))
    assemblieName = list(dummyData)[0]
    assemblieLinks = dummyData[assemblieName]
    platypus = ncbiData(assemblieName, assemblieLinks)

if __name__ == "__main__":
    test_ncbiDataInit()
