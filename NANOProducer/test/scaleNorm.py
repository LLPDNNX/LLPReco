import os
import sys
import uproot
import re
import json
import numpy as np
import math

scaleNorm = {}


for year in ['2016','2017','2018']:
    basepath = "/vols/cms/mkomm/HNL/CMSSW_10_2_22/src/files_201117/"+year

    txtFiles = sorted(os.listdir(basepath))



    for txtFile in txtFiles:
        process = txtFile.split(".")[0]
        
        
        if not (
            re.match("TT\S+powheg\S+",process)>=0 or \
            process.find("WJetsToLNu")>=0 or \
            re.match("W\S+amcatnlo\S+",process)>=0 or \
            process.find("TTToHadronic")>=0 or \
            process.find("TTToSemiLeptonic")>=0 or \
            process.find("TTTo2L2Nu")>=0 or \
            process.find("ST_")>=0 or \
            re.match("DYJets\S+amcatnlo\S+",process) or \
            re.match("DY\S+",process) or \
            re.match("WG\S+",process)>=0 or \
            re.match("ZG\S+",process)>=0 
        ):
            print "skip ",process
            continue
            
        sumNominal = 0.
        sumUp = 0.
        sumDown = 0.
            
        f = open(os.path.join(basepath,txtFile))
        print "reading ... ",txtFile
        for l in f:
            if len(l)>0:
                fileName = l.replace("\n","").replace("\r","")
                
                f = uproot.open(fileName)
                data = f["Events"].arrays(["LHEScaleWeight","nLHEScaleWeight"])
                if data["nLHEScaleWeight"][0]<9:
                    print "Process %s has no scale weights stored"%(process)
                    sumNominal = -1
                    break
                
                #NB: weights are always positive
                for iev in range(len(data['nLHEScaleWeight'])):
                    
                    nominal = data['LHEScaleWeight'][iev][4]
                    up = data['LHEScaleWeight'][iev][8]
                    down = data['LHEScaleWeight'][iev][0]
                    if iev<10:
                        print "%.3f, %.3f, %.3f"%(nominal,up,down)
                    if math.fabs(nominal-1)>1 or math.fabs(up-1)>1 or math.fabs(down-1)>1:
                        nominal = 1
                        up = 1
                        down = 1
                        
                    sumNominal += nominal
                    sumDown += down
                    sumUp += up
                break
                if sumNominal>1e5:
                    break
                
        if sumNominal>0:
            print "%.3f, %.3f"%(sumUp/sumNominal,sumDown/sumNominal)
            scaleNorm[process] = {"up": sumUp/sumNominal, "down": sumDown/sumNominal}
        
with open("scaleNormTable.json","w") as f:
    json.dump(
        scaleNorm,
        f,
        ensure_ascii=True, 
        check_circular=True, 
        allow_nan=True, 
        cls=None, 
        indent=4, 
        sort_keys=True, 
    )
        

        
