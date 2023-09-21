import sys
from ncbiData import ncbiData

if __name__ == "__main__":
    species = ncbiData(sys.argv[1])
    species.runPremrnaBlast()
    species.runMrnaBlast()
