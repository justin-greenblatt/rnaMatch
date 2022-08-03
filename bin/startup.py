from subprocess import Popen, PIPE
from settings import sConfig, dConfig, pConfig
import os


#formatDisk
formatDiskCommand = ["sudo", "mkfs.ext4", "-m", "0", "-E", "lazy_itable_init=0,lazy_journal_init=0,discard", "/dev/sdb"]

pFormat = Popen(formatDiskCommand, stdout = PIPE, stdin = PIPE)
pFormat.wait()
print("-----formatedDisk-----")

def createDir(dirName):
    if not os.path.isdir(dirName):
        os.mkdir(dirName)

createDir(dConfig["common"]["data_path"])
mountDiskCommand = ["sudo", "mount", "-o", "discard,defaults", "/dev/sdb", dConfig["common"]["data_path"]]


pMount = Popen(mountDiskCommand, stdout = PIPE, stdin = PIPE)
pMount.wait()
print("-----mounted disk-----")

for d in dConfig["resources"].values():
    createDir(d)
    print("---created directorie : "+d)
