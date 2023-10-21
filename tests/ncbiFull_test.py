import sys, os, json
sys.path.insert(1, os.path.join(os.environ.get("HOME"), "rnaMatch/bin"))
from ncbiData import ncbiData
import logging
import sys
class TestNcbiData:
    print(f"Extracting from {sys.argv[1]}")
    i = ncbiData(sys.argv[1])

    def test_full(self):
        TestNcbiData.i.getResource("genome")
        TestNcbiData.i.getResource("gtf")
        TestNcbiData.i.getResource("mrna")
        TestNcbiData.i.runMrnaBlast()
        TestNcbiData.i.runPremrnaBlast()
        TestNcbiData.i.compressAndUpload(['mrna',
                                 'mrna_blast_control_out_summary',
                                 'premrna_blast_test_out',
                                 'mrna_blast_control_out',
                                 'mrna_blast_test_out_summary',
                                 'premrna_blast_test_out_summary',
                                 'premrna_blast_control_out_summary',
                                 'premrna',
                                 'premrna_blast_control_out'])

        TestNcbiData.i.upload()

if __name__ == "__main__":
    T = TestNcbiData()
    T.test_full()