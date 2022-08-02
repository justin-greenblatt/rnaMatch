from subprocess import Popen, PIPE
import sys, os
sys.path.insert(1, os.path.join(os.environ.get("HOME"), "blastWeb/bin"))
from settings.directories import PROCESSES_CONF, REV_BLAST_PATH, TEST_REV_BLAST_IN_FILE, TEST_REV_BLAST_CONTROL_OUT_FILE, TEST_REV_BLAST_OUT_FILE
from configparser import ConfigParser

pConfig = ConfigParser()
pConfig.read(PROCESSES_CONF)

class TestRevBlast:
        
    command  = list([a for a in [
                    pConfig["revBlast"]["user"],
                    pConfig["revBlast"]["userFlags"],
                    pConfig["revBlast"]["interpreter"],
                    pConfig["revBlast"]["interpreterFlags"]]
                    if a])

    def test_revBlatsPlus(self):
        testCommand = TestRevBlast.command + [REV_BLAST_PATH, TEST_REV_BLAST_IN_FILE, TEST_REV_BLAST_CONTROL_OUT_FILE, "plus"]
        print(testCommand)
        p = Popen(testCommand, stdout = PIPE)
        assert p.wait() == 0

    def test_revBlastMinus(self):
        testCommand = TestRevBlast.command + [REV_BLAST_PATH, TEST_REV_BLAST_IN_FILE, TEST_REV_BLAST_OUT_FILE, "minus"]
        print(testCommand)
        p = Popen(testCommand, stdout = PIPE)
        assert p.wait() == 0

