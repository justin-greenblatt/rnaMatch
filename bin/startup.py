from subprocess import Popen, PIPE
import os
from settings import dConfig as homelessConfig
from settings import configsPath

#configure home Directorie
homelessConfig["common"]["HOME_DIR"] = os.environ.get("HOME")
dConfigOut = open('/' + str(os.path.join(*os.path.realpath(__file__).split('/')[:-1], "settings", "directories.ini")), 'w')
homelessConfig.write(dConfigOut)
dConfigOut.close()

from settings import sConfig, dConfig, pConfig

def runCommand(command):
    p = Popen(command, stdout = PIPE, stdin = PIPE)
    print(" ".join(command))
    p.wait()
    print("finished")
    return p.communicate()

for p in sConfig["apt"].values():
    runCommand(["sudo", "apt", "install", p, "-y"])
for q in sConfig["pip"].values():
    runCommand(["sudo", "pip3", "install", q, "--no-input"])

def createDir(dirName):
    if not os.path.isdir(dirName):
        os.mkdir(dirName)       
        runCommand(["chmod", "-R", "777", dirName])

createDir(dConfig["common"]["out_path"])
createDir(dConfig["common"]["data_path"])
mountDiskCommand = ["sudo", "mount", sConfig["nfs"]["IP"], dConfig["common"]["nfs_path"]]
c = runCommand(mountDiskCommand)
print(c)
for f in dConfig["resourceFolders"].values():
    createDir(f)

for d in dConfig["resources"].values():
    createDir(d)

