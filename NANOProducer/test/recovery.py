import CRABClient
from CRABClient.UserUtilities import config, getLumiListInValidFiles
from WMCore.DataStructs.LumiList import LumiList
import json, urllib

failedTasks = {
        "SingleMuon_Run2017F": "2017",
        #"EGamma_Run2018A": "2018",
        #"SingleElectron_Run2017C": "2017",
        #"SingleMuon_Run2017E": "2017"
        }

taskDAS = {
        "SingleMuon_Run2017F": "/SingleMuon/vcepaiti-SingleMuon_Run2017F-f34f4fd09de6de3dea965fd010802907/USER"
        #"EGamma_Run2018A": "/EGamma/vcepaiti-EGamma_Run2018A-2c4386821e9589ec1e50325b77b4b199/USER",
        #"SingleElectron_Run2017C": "/SingleElectron/vcepaiti-SingleElectron_Run2017C-f34f4fd09de6de3dea965fd010802907/USER",
        #"SingleMuon_Run2017E": "/SingleMuon/vcepaiti-SingleMuon_Run2017E-f34f4fd09de6de3dea965fd010802907/USER"
        }
for task, year in failedTasks.items():
    if year == '2016':
        url='https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification//Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
    if year == '2017':
        url='https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
    if year == '2018' or year== '2018D':
        url='https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt' 

    json_url = urllib.urlopen(url)
    officialLumiMaskJson = json.loads(json_url.read())
    print(officialLumiMaskJson)
    officialLumiMask = LumiList(compactList=officialLumiMaskJson)

    failedTaskLumis = getLumiListInValidFiles(dataset=taskDAS[task], dbsurl='phys03')
    newLumiMask = officialLumiMask - failedTaskLumis

    newLumiMask.writeJSON('{}.json'.format(task))


