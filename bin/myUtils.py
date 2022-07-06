"""
Developed for python3.9
justingreenblatt@github.com | 09/05/2022
Utils used in more than one object.
"""
import requests
from subprocess import Popen, PIPE
from os import listdir, getcwd
from shutil import unpack_archive
from os import listdir
import logging
from time import time
from settings.logs import DOWNLOADS_LOG_PATH, DOWNLOADS_LOG_LEVEL
from typing import List, Dict, TypedDict
import re

def downloadFromURL(url : str, filename : str = False, decompress : bool = False) -> str:
    """
    Download data from a url link to a file.
    if decompress = True, then use gunzip on the file.
    IMPORTANT!!! it returns the file path of downloaded and processed content
    Logging is in  downloads
    """

    #dissable warnings
    requests.packages.urllib3.disable_warnings(requests.packages.urllib3.exceptions.InsecureRequestWarning)
    
    if not filename:
        filename = url.split('/')[-1]
    h = open(filename, 'wb+')

    r = requests.get(url, stream=True, allow_redirects=True, verify=False)
    data = r.content
    h = open(filename, 'wb+')
    h.write(data)
    h.close()

    if decompress:
        #MODIFY LATER here for modular decompression
        command = ["gzip","-v", "-d", filename]
        dP = Popen(command, stdin = PIPE, stdout = PIPE)
        dP.wait()
        filename = filename.rstrip(".gz")
    return filename

def getLinks(link : str) -> List[str]:
    return re.findall(r'href=\"(.*?)\"', requests.get(link, verify = False).text)

class dictCompareOut(TypedDict):
            added: Dict[str, str]
            removed: Dict[str, str]
            altered: Dict[str, TypedDict('Changed', {'old': str, 'new': str})]
            same: Dict[str,str]

def dictComparison(old : Dict[str,str], new: Dict[str,str]) -> dictCompareOut:
    out = {}
    out["added"] = {k : new[k] for k in list(set(new.keys()) - set(old.keys()))} 
    out["removed"] = {k: old[k] for k in list(set(old.keys()) - set(new.keys()))}
    out["altered"] = {k : {"old" : old[k], "new" : new[k]} for k in list(set(old.keys()) & set(new.keys())) if old[k] != new[k]}
    out["same"] = {k : old[k] for k in list(set(old.keys()) & set(new.keys())) if old[k] == new[k]}
    return out
