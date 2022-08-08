import sys, os

sys.path.insert(0, os.path.join(os.environ.get("HOME"), "blastWeb/bin"))
from subprocess import Popen, PIPE
from settings import dConfig, pConfig

def test_genomeWalk():
    if os.path.isfile(dConfig["tests"]["TEST_GENOME_WALK_OUT"]):
        os.remove(dConfig["tests"]["TEST_GENOME_WALK_OUT"])
    if os.path.isfile(dConfig["tests"]["TEST_GENOME_WALK_CONTROL_OUT"]):
        os.remove(dConfig["tests"]["TEST_GENOME_WALK_CONTROL_OUT"])


    command  = list([a for a in [
            pConfig["genomeWalk"]["USER"],
            pConfig["genomeWalk"]["USER_FLAGS"],
            pConfig["genomeWalk"]["INTERPRETER"],
            pConfig["genomeWalk"]["INTERPRETER_FLAGS"],
            dConfig["scripts"]["GENOME_WALK_PATH"], 
            dConfig["tests"]["TEST_GENOME"], 
            dConfig["tests"]["TEST_GTF"],
            dConfig["tests"]["TEST_GENOME_WALK_OUT"],
            dConfig["tests"]["TEST_GENOME_WALK_CONTROL_OUT"]]
                        if a and a != "None"])
    print(command)
    p = Popen(command, stdout = PIPE)
    assert p.wait() == 0
