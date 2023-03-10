import subprocess
import re
import time
import os

from fcntl import fcntl, F_GETFL, F_SETFL
from os import O_NONBLOCK, read

#basepath = "srm://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mkomm/HNL/NANOX_201117"
#basepath = "srm://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/vcepaiti/HNL/NANOX_201117"
basepath = "srm://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mkomm/HNL/NANOX_221221"

def listdir(path):
    #print "query path ... ",path
    ret = -1
    retry = 0
    fetchedOutput = ""
    while (ret!=0 and retry<20):
        if retry>0:
            print "retry ",retry,"/ 20 ...",path
            time.sleep(1+2*retry)
        proc = subprocess.Popen([
                "gfal-ls",
                "-t",
                "90", #increase default timeout
                "-l",
                path
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        flags = fcntl(proc.stdout, F_GETFL)
        fcntl(proc.stdout, F_SETFL, flags | O_NONBLOCK)
        waiting = 200+10*retry
        while (proc.poll()==None and waiting>0):
            time.sleep(0.1)
            waiting-=1
            try:
                fetchedOutput+=read(proc.stdout.fileno(), 1024)
            except OSError:
                pass
        if (proc.poll()==None):
            proc.terminate()
            retry+=1
            continue
        ret = proc.returncode
        if ret!=0:
            linesOut,linesErr = proc.communicate()
            print linesErr
        retry+=1
        
    linesOut,linesErr = proc.communicate()
    if ret==0:
        linesOut = filter(None,(fetchedOutput+linesOut).replace('\t','').split('\n'))
        fileList = []
        for l in linesOut:
            
            items = filter(None,l.split(" "))
            if len(items)!=9:
                print "Error while parsing items: ",items
            fileList.append({
                "isDir":items[0].startswith('d'),
                "permission":items[0],
                "nlinks":items[1],
                "group":items[2],
                "user":items[3],
                "size":items[4]+" "+items[5],
                "date":items[6],
                "time":items[7],
                "name":items[8]
            })
        return fileList
    else:
        print "An error occured: ",linesErr
        return []
        
def listAll(path,pattern,depth=10):
    #print "query path ... ",path
    proc = subprocess.Popen([
            "srmls",
            "-recursion_depth="+str(depth),
            path
        ],
        stdout=subprocess.PIPE
    )
    proc.wait()
    ret = proc.returncode
    linesOut,linesErr = proc.communicate()
    if ret==0:
        linesOut = filter(None,linesOut.replace('\t','').split('\n'))
        fileList = []
        for l in linesOut:
            
            items = filter(None,l.split(" "))
            if len(items)!=2:
                print "Error while parsing items: ",items
            if pattern.match(items[1]):
                fileList.append(items[1])
        return fileList
    else:
        print "An error occured: ",linesErr
        return []
        
def findFiles(path,pattern):
    matches = []
    for f in listdir(path):
        if f["isDir"]:
            if f["name"].find("log")>=0 or f["name"].find("failed")>=0:
                continue
            matches.extend(findFiles(path+"/"+f["name"],pattern))
        else:
            if pattern.match(path+"/"+f["name"]):
                matches.append(path+"/"+f["name"])
    return matches
    
folders = listdir(basepath)


for i,folder in enumerate(folders):
    if not folder["isDir"]:
        continue
    '''   
    if folder["name"].find("Single")<0 and folder["name"].find("EG")<0:
        continue
    '''
    '''
    if folder["name"].find("ST_")<0:
        continue
    '''
    print i+1,"/",len(folders),folder["name"]
    
    fileList = findFiles(
        basepath+'/'+folder["name"],
        re.compile("\S+_[0-9]+.root")
    )
    if len(fileList)==0:
        print "  WARNING: directory empty -> skip"
        continue
        
    filesPerDate = {}
    for f in fileList:
        g = re.match("\\S+([0-9][0-9][0-9][0-9][0-9][0-9])_([0-9][0-9][0-9][0-9][0-9][0-9])/[0-9]+/\\S+_[0-9]+.root",f)
        if not g:
            raise Exception("File path not recognized: "+f)
        #print g.groups()[0],g.groups()[1],f
        k = g.groups()[0]+'_'+g.groups()[1]
        if not filesPerDate.has_key(k):
            filesPerDate[k] = []
        filesPerDate[k].append(f)
        
    if len(filesPerDate.keys())>1:
        print "  WARNING - found multiple jobs: "
        for k in filesPerDate.keys():
            print '    ',k,len(filesPerDate[k])
    kLatest = sorted(filesPerDate.keys(),key=lambda x: int(x[0:6])*1e8+int(x[7:]))[-1]
    
    #print kLatest
    fileList = filesPerDate[kLatest]
    
    print '   #files=',len(fileList)
    
    outFile = open(folder["name"]+".txt","w")
    for f in fileList:
        if f.find("failed")>=0:
            continue
        #print f
        outFile.write(f.replace("srm","root")+"\n")
    outFile.close()
   
    
