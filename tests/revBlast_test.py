from subprocess import Popen, PIPE
import sys, os

sys.path.insert(0, os.path.join(os.environ.get("HOME"), "blastWeb/bin"))

from settings import pConfig, dConfig

class TestRevBlast:

    if os.path.isfile(dConfig["tests"]["TEST_REV_BLAST_CONTROL_OUT_FILE"]):
        os.remove(dConfig["tests"]["TEST_REV_BLAST_CONTROL_OUT_FILE"])

    command  = list([a for a in [
                    pConfig["revBlast"]["USER"],
                    pConfig["revBlast"]["USER_FLAGS"],
                    pConfig["revBlast"]["INTERPRETER"],
                    pConfig["revBlast"]["INTERPRETER_FLAGS"]]
                    if a])

    def test_revBlatsPlus(self):
        testCommand = TestRevBlast.command + [dConfig["scripts"]["REV_BLAST_PATH"],
                                              dConfig["tests"]["TEST_REV_BLAST_IN_FILE"],
                                              dConfig["tests"]["TEST_REV_BLAST_CONTROL_OUT_FILE"],
                                              "plus"]
        print(testCommand)
        p = Popen(testCommand, stdout = PIPE)
        assert p.wait() == 0

    def test_revBlastMinus(self):

        if os.path.isfile(dConfig["tests"]["TEST_REV_BLAST_OUT_FILE"]):
            os.remove(dConfig["tests"]["TEST_REV_BLAST_OUT_FILE"])

        testCommand = TestRevBlast.command + [dConfig["scripts"]["REV_BLAST_PATH"],
                dConfig["tests"]["TEST_REV_BLAST_IN_FILE"],
                dConfig["tests"]["TEST_REV_BLAST_OUT_FILE"],
                "minus"]
        p = Popen(testCommand, stdout = PIPE)
        assert p.wait() == 0
