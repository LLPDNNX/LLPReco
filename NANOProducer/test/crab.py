from WMCore.Configuration import Configuration
import datetime,sys,os
import copy
import math
import urllib, json
from CRABClient.UserUtilities import getUsernameFromSiteDB, getLumiListInValidFiles


myJobs = {
    "TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8-2016":{
        "inputDataset":"/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    "SingleMuon_Run2016D":{
        "inputDataset": "/SingleMuon/Run2016D-17Jul2018-v1/MINIAOD",
        "isData": True,
        "year": 2016
    },
    "GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8":{
        "inputDataset": "/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/kjpena-RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_MINIAODSIM-bd3e7bcff6c9bcad356ea4ed7e4f08b4/USER",
        "year": 2016,
        "LLPtype": "hss"
    }
}


'''
    "TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-2016":{
        "inputDataset":"/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    "TTJets_TuneCP5_13TeV-madgraphMLM-pythia8-2017":{
        "inputDataset":"/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    "TTJets_TuneCP5_13TeV-madgraphMLM-pythia8-2018":{
        "inputDataset":"/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
    "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_ext1-2017":{
        "inputDataset":"/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2

    },
    "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8-2018":{
        "inputDataset":"/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },

    "TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8-2016":{
        "inputDataset":"/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    "TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8-2017":{
        "inputDataset":"/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8-2018":{
        "inputDataset":"/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
'''

'''
myJobs = {
    "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-2016":{
        "inputDataset":"/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v3/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-ext1-2016":{
        "inputDataset":"/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-ext2-2016":{
        "inputDataset":"/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    
    
    "DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2017":{
        "inputDataset":"/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    "DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2017":{
        "inputDataset":"/DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    "DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2017":{
        "inputDataset":"/DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    
    
    "DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018":{
        "inputDataset":"/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
    "DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018":{
        "inputDataset":"/DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
    "DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018":{
        "inputDataset":"/DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
}
'''



'''
myJobs = {

    "WToLNu_0J_13TeV-amcatnloFXFX-pythia8-2016":{
        "inputDataset":"/WToLNu_0J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_backup_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    "WToLNu_0J_13TeV-amcatnloFXFX-pythia8-ext1-2016":{
        "inputDataset":"/WToLNu_0J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2017":{
        "inputDataset":"/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    
    "WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018":{
        "inputDataset":"/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
    
    
    
    "WToLNu_1J_13TeV-amcatnloFXFX-pythia8-2016":{
        "inputDataset":"/WToLNu_1J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_backup_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8-ext1-2017":{
        "inputDataset":"/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    
    "WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018":{
        "inputDataset":"/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
    
    
    
    "WToLNu_2J_13TeV-amcatnloFXFX-pythia8-ext4-2016":{
        "inputDataset":"/WToLNu_2J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext4-v1/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    
    "WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2017":{
        "inputDataset":"/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },

    "WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018":{
        "inputDataset":"/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
}
'''
'''
myJobs = {
    "QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    "QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    "QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    "QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia-2016":{
        "inputDataset":"/QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },




    "QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    "QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    "QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    "QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    "QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    "QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    "QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    "QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    "QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    "QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    "QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },
    
    
    
    "QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
    "QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v4/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
    "QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
    "QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
    "QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
    "QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
    "QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
    "QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext3-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
    "QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
    "QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
    "QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext3-v2/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
    "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },
}
'''



'''
myJobs = {}
for folder in os.listdir('hnl_mu'):
    name = folder.replace('.txt','')
    name = name.replace(".","_") #'.' not allowed in name
    name = name.replace("massiveAndCKM_LO_","")
    if name.find('_Dirac_')<0:
        continue
    if name.find('_pre2017_')>=0:
        continue
    if name.find('Moriond17_aug2018')>0 and name.find('miniAODv3')<0:
        continue
        
    print name
    
    year = None
    if name.find('Moriond17_aug2018')>0:
        year = 2016
    elif name.find('Fall17')>0:
        year = 2017
    elif name.find('Autumn18')>0:
        year = 2018
    if not year:
        print "Year not found for ",folder
        continue
    name = name+"-"+str(year)
        
    userInputFiles = os.path.abspath(os.path.join('hnl_mu',folder))
        
    if len(open(userInputFiles).readlines())<2:
        print "WARNING: only 1 file found"
        continue
    
    myJobs[name] = {
        "userInputFiles": userInputFiles,
        "year": year,
        "unitsPerJob":20,
        "isData":False,
        'whitelist': ['T2_BE_IIHE','T2_BE_UCL','T2_CH_CERN','T2_UK_London_IC'],
    }
    #break

'''

'''
myJobs = {
    "SingleMuon_Run2016B_ver1":{
        "inputDataset": "/SingleMuon/Run2016B-17Jul2018_ver1-v1/MINIAOD",
        "isData": True,
        "year": '2016'
    },
    "SingleMuon_Run2016B_ver2":{
        "inputDataset": "/SingleMuon/Run2016B-17Jul2018_ver2-v1/MINIAOD",
        "isData": True,
        "year": '2016'
    },
    "SingleMuon_Run2016C":{
        "inputDataset": "/SingleMuon/Run2016C-17Jul2018-v1/MINIAOD",
        "isData": True,
        "year": '2016'
    },
    "SingleMuon_Run2016D":{
        "inputDataset": "/SingleMuon/Run2016D-17Jul2018-v1/MINIAOD",
        "isData": True,
        "year": '2016'
    },
    "SingleMuon_Run2016E":{
        "inputDataset": "/SingleMuon/Run2016E-17Jul2018-v1/MINIAOD",
        "isData": True,
        "year": '2016'
    },
    "SingleMuon_Run2016F":{
        "inputDataset": "/SingleMuon/Run2016F-17Jul2018-v1/MINIAOD",
        "isData": True,
        "year": '2016'
    },
    "SingleMuon_Run2016G":{
        "inputDataset": "/SingleMuon/Run2016G-17Jul2018-v1/MINIAOD",
        "isData": True,
        "year": '2016'
    },
    "SingleMuon_Run2016H":{
        "inputDataset": "/SingleMuon/Run2016H-17Jul2018-v1/MINIAOD",
        "isData": True,
        "year": 2016
    },
    
    
    
    "SingleMuon_Run2017B":{
        "inputDataset": "/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD",
        "isData": True,
        "year": '2017'
    },
    "SingleMuon_Run2017C":{
        "inputDataset": "/SingleMuon/Run2017C-31Mar2018-v1/MINIAOD",
        "isData": True,
        "year": 2017
    },
    "SingleMuon_Run2017D":{
        "inputDataset": "/SingleMuon/Run2017D-31Mar2018-v1/MINIAOD",
        "isData": True,
        "year": '2017'
    },
    "SingleMuon_Run2017E":{
        "inputDataset": "/SingleMuon/Run2017E-31Mar2018-v1/MINIAOD",
        "isData": True,
        "year": '2017'
    },
    "SingleMuon_Run2017F":{
        "inputDataset": "/SingleMuon/Run2017F-31Mar2018-v1/MINIAOD",
        "isData": True,
        "year": '2017'
    },
    
    
    "SingleMuon_Run2018A":{
        "inputDataset": "/SingleMuon/Run2018A-17Sep2018-v2/MINIAOD",
        "isData": True,
        "year": '2018'
    },
    "SingleMuon_Run2018B":{
        "inputDataset": "/SingleMuon/Run2018B-17Sep2018-v1/MINIAOD",
        "isData": True,
        "year": '2018'
    },
    "SingleMuon_Run2018C":{
        "inputDataset": "/SingleMuon/Run2018C-17Sep2018-v1/MINIAOD",
        "isData": True,
        "year": '2018'
    },
    
    
    "SingleMuon_Run2018D":{
        "inputDataset": "/SingleMuon/Run2018D-22Jan2019-v2/MINIAOD",
        "isData": True,
        "year": '2018D'
    }
}
'''

'''
hss2016Files = [
'/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-0_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-0p05_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-0p1_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-10_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-1_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-2000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-25_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-5000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-500_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-50_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-5_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-20_ctauS-0_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-20_ctauS-0p05_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-20_ctauS-0p1_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-20_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-20_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-20_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-20_ctauS-10_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-20_ctauS-1_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-20_ctauS-2000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-20_ctauS-25_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-20_ctauS-5000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-20_ctauS-500_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-20_ctauS-50_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-20_ctauS-5_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-25_ctauS-0_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-25_ctauS-0p05_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-25_ctauS-0p1_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-25_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-25_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-25_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-25_ctauS-10_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-25_ctauS-1_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-25_ctauS-25_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-25_ctauS-5000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-25_ctauS-500_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-25_ctauS-50_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-25_ctauS-5_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-30_ctauS-0_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-30_ctauS-0p05_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-30_ctauS-0p1_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-30_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-30_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-30_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-30_ctauS-10_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-30_ctauS-1_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-30_ctauS-2000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-30_ctauS-25_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-30_ctauS-5000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-30_ctauS-500_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-30_ctauS-50_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-30_ctauS-5_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-40_ctauS-0_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-40_ctauS-0p05_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-40_ctauS-0p1_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-40_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-40_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-40_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-40_ctauS-10_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-40_ctauS-1_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-40_ctauS-2000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-40_ctauS-25_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-40_ctauS-5000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-40_ctauS-500_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-40_ctauS-50_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-40_ctauS-5_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-50_ctauS-0_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-50_ctauS-0p05_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-50_ctauS-0p1_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-50_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-50_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-50_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-50_ctauS-10_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-50_ctauS-1_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-50_ctauS-2000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-50_ctauS-25_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-50_ctauS-5000_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-50_ctauS-500_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-50_ctauS-50_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
'/GluGluH_HToSSTobbbb_MH-125_MS-50_ctauS-5_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets-1094e4be9b2393641646b907f4de5a1a/USER',
]

myJobs = {}

for hss2016File in hss2016Files:
    name = hss2016File.split('/')[1]
    myJobs[name] = {
        "inputDataset":hss2016File,
        "year": 2016,
        "unitsPerJob": 2,
        'inputDBS':'prod/phys03'
    }
'''

requestName = "NANOX_060320"
userName = getUsernameFromSiteDB() 
configTmpl = Configuration()

configTmpl.section_('General')
configTmpl.General.transferOutputs = True
configTmpl.General.transferLogs = False

configTmpl.section_('JobType')
configTmpl.JobType.psetName = "LLPReco/NANOProducer/test/produceNANO.py"
configTmpl.JobType.pluginName = 'Analysis'
configTmpl.JobType.outputFiles = ['nano.root']
configTmpl.JobType.allowUndistributedCMSSW = True
configTmpl.JobType.maxJobRuntimeMin= 20*60
configTmpl.JobType.pyCfgParams = []
configTmpl.JobType.inputFiles = []
configTmpl.JobType.maxMemoryMB = 2499
configTmpl.section_('Data')
configTmpl.Data.useParent = False
configTmpl.section_('Site')
configTmpl.Site.storageSite = 'T2_UK_London_IC'

'''
configTmpl.Data.ignoreLocality = True
configTmpl.Site.whitelist = [
    'T2_BE_IIHE','T2_BE_UCL','T2_CH_CERN','T2_UK_London_IC',
    'T2_DE_DESY','T2_DE_RWTH',
]
'''

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    from multiprocessing import Process

    def submit(config):
        try:
            crabCommand('submit',  config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################


    for i,jobName in enumerate(sorted(myJobs.keys())):

        isData = False
        myJob = myJobs[jobName]
        i=i+1
        config = copy.deepcopy(configTmpl)
        config.General.requestName = jobName+"_"+requestName
        print config.General.requestName
        config.General.workArea = "crab/"+requestName+"/"+jobName
        config.Data.outLFNDirBase = "/store/user/"+userName+"/LLP/"+requestName+"/"+jobName
        userInputFiles = myJob.get('userInputFiles', None)
        if not userInputFiles:
            config.Data.inputDataset = myJob["inputDataset"]
            config.Data.inputDBS = myJob.get('inputDBS', 'global')
            config.Data.publication = True
        else:
            config.Data.userInputFiles = map(lambda f: f.replace('\n','').replace('\r',''),open(userInputFiles).readlines())
            config.Data.outputPrimaryDataset = jobName
            config.Data.publication = True

        config.Data.outputDatasetTag = jobName
            
        isData = myJob.get('isData', False)
        year = str(myJob.get('year', '0'))
        
        if year not in ['2016','2017','2018','2018D']:
            print "ERROR: Year invalid: ",year
            continue
        print "year:", year
        
        config.JobType.pyCfgParams.append("year="+str(year))

        if isData:
            config.JobType.pyCfgParams.append("isData=True")
            config.Data.splitting = 'LumiBased'
            config.Data.unitsPerJob = myJob.get('unitsPerJob', 80)
            if year == '2016':
                config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification//Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
            if year == '2017':
                config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
            if year == '2018' or year== '2018D':
                config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt' 
        else:
            config.JobType.pyCfgParams.append("isData=False")
            config.Data.splitting = 'FileBased'
            config.Data.unitsPerJob = myJob.get('unitsPerJob', 1)

        if "params" in myJob:
            params = myJob["params"]

            for param in params:
                config.JobType.pyCfgParams.append(str(param))
                
        if "whitelist" in myJob:
            config.Site.whitelist = myJob['whitelist']
            
        if "blacklist" in myJob:
            config.Site.whitelist = myJob['blacklist']

        
        
        if not os.path.exists(configTmpl.JobType.psetName):
            print "\nConfiguration file ", pSet, "does not exist.  Aborting..."
            sys.exit(1)
        
        
        
        if os.path.isdir(os.path.join(os.getcwd(),config.General.workArea)):
            print "Output directory ",os.path.join(os.getcwd(),config.General.workArea)," exists -> skipping"
            print
            continue
            
        print config,
            
        print "Submitting job ",i," of ",len(myJobs.keys()),":",config.General.workArea

        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        print
        print
        
