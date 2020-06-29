import CRABClient
import datetime,sys,os
import copy
import math
import urllib, json
from WMCore.Configuration import Configuration
#from CRABClient.UserUtilities import getUsernameFromSiteDB, getLumiListInValidFiles

myJobs = {
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

    "TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8-2016":{
        "inputDataset":"/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8-2016":{
        "inputDataset":"/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },


    "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_ext1-2017":{
        "inputDataset":"/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },

    "TTToHadronic_TuneCP5_13TeV-powheg-pythia8-2017":{
        "inputDataset":"/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },

    "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8-2018":{
        "inputDataset":"/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },

    "TTToHadronic_TuneCP5_13TeV-powheg-pythia8-2018":{
        "inputDataset":"/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
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

    
    "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1-2016":{
        "inputDataset":"/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM",
        "year": 2016
    },

    "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1-2016":{
        "inputDataset":"/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM",
        "year": 2016
    },

    "ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8-2017":{
        "inputDataset":"/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017
    },
    
    "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8-2017":{
        "inputDataset":"/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
        "year": 2017
    },

    "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8-2018":{
        "inputDataset":"/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM",
        "year": 2018
    },

    "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8-2018":{
        "inputDataset":"/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM",
        "year": 2018
    },
    
    "ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1-2016":{
        "inputDataset":"/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
        "year": 2016
    },
    
    "ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8-2016":{
        "inputDataset":"/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
        "year": 2016
    },
       
    "ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8-2017":{
        "inputDataset":"/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017
    },
    
    "ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8-2017":{
        "inputDataset":"/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017
    },
    
    "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8-2018":{
        "inputDataset":"/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018
    },
    
    "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8-2018":{
        "inputDataset":"/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018
    },

    "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-2016":{
        "inputDataset":"/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-2016":{
        "inputDataset":"/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v3/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-ext2-2016":{
        "inputDataset":"/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8-2017":{
        "inputDataset":"/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 3
    },

    "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8-2017":{
        "inputDataset":"/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 3
    },

    "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8-2018":{
        "inputDataset":"/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 3
    },

    "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018": {
        "inputDataset":"/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 3
    },

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

    "ZGToLLG_01J_5f_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-2016":{
        "inputDataset":"/ZGToLLG_01J_5f_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8-2017":{
        "inputDataset": "/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "ZGToLLG_01J_5f_PDFWeights_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018":{
        "inputDataset": "/ZGToLLG_01J_5f_PDFWeights_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-2016":{
        "inputDataset": "/WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8-2017":{
        "inputDataset": "/WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },

    "WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8-2018":{
        "inputDataset": "/WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },

    "WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8-2016":{
        "inputDataset": "/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8-2017":{
        "inputDataset": "/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },


    "WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018":{
        "inputDataset": "/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },


    "TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8-2016":{
        "inputDataset": "/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8-2017":{
        "inputDataset": "/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },


    "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8-2018":{
        "inputDataset": "/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },


    "TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8-2016":{
        "inputDataset": "/TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8-2017":{
        "inputDataset": "/TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },

    "TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8-2018":{
        "inputDataset": "/TGJets_TuneCP5_13TeV_amcatnlo_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },

    "QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },
    "QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob": 2
    },

    "QCD_Pt-15to20_EMEnriched_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-15to20_EMEnriched_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },

    "QCD_Pt-20to30_EMEnriched_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-20to30_EMEnriched_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },


    "QCD_Pt-30to50_EMEnriched_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-30to50_EMEnriched_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },


    "QCD_Pt-50to80_EMEnriched_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-50to80_EMEnriched_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },


    "QCD_Pt-80to120_EMEnriched_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-80to120_EMEnriched_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },


    "QCD_Pt-120to170_EMEnriched_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-120to170_EMEnriched_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },


    "QCD_Pt-170to300_EMEnriched_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-170to300_EMEnriched_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },

    "QCD_Pt-300toInf_EMEnriched_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt-300toInf_EMEnriched_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 2
    },


    "QCD_Pt-15to20_EMEnriched_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-15to20_EMEnriched_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },

    "QCD_Pt-20to30_EMEnriched_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-20to30_EMEnriched_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },

    "QCD_Pt-30to50_EMEnriched_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-30to50_EMEnriched_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },

    "QCD_Pt-50to80_EMEnriched_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-50to80_EMEnriched_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },

    "QCD_Pt-80to120_EMEnriched_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-80to120_EMEnriched_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },

    "QCD_Pt-120to170_EMEnriched_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-120to170_EMEnriched_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },

    "QCD_Pt-170to300_EMEnriched_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-170to300_EMEnriched_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },


    "QCD_Pt-300toInf_EMEnriched_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt-300toInf_EMEnriched_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 2
    },


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
    "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016":{
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
myJobs = {
    "SingleElectron_Run2016B_ver2":{
        "inputDataset": "/SingleElectron/Run2016B-17Jul2018_ver2-v1/MINIAOD",
        "isData": True,
        "year": '2016'
    },
    "SingleElectron_Run2016C":{
        "inputDataset": "/SingleElectron/Run2016C-17Jul2018-v1/MINIAOD",
        "isData": True,
        "year": '2016'
    },
    "SingleElectron_Run2016D":{
        "inputDataset": "/SingleElectron/Run2016D-17Jul2018-v1/MINIAOD",
        "isData": True,
        "year": '2016'
    },
    "SingleElectron_Run2016E":{
        "inputDataset": "/SingleElectron/Run2016E-17Jul2018-v1/MINIAOD",
        "isData": True,
        "year": '2016'
    },
    "SingleElectron_Run2016F":{
        "inputDataset": "/SingleElectron/Run2016F-17Jul2018-v1/MINIAOD",
        "isData": True,
        "year": '2016'
    },
    "SingleElectron_Run2016G":{
        "inputDataset": "/SingleElectron/Run2016G-17Jul2018-v1/MINIAOD",
        "isData": True,
        "year": '2016'
    },
    "SingleElectron_Run2016H":{
        "inputDataset": "/SingleElectron/Run2016H-17Jul2018-v1/MINIAOD",
        "isData": True,
        "year": '2016'
    },
        
    "SingleElectron_Run2017B":{
        "inputDataset": "/SingleElectron/Run2017B-31Mar2018-v1/MINIAOD",
        "isData": True,
        "year": '2017'
    },
    "SingleElectron_Run2017C":{
        "inputDataset": "/SingleElectron/Run2017C-31Mar2018-v1/MINIAOD",
        "isData": True,
        "year": '2017'
    },
    "SingleElectron_Run2017D":{
        "inputDataset": "/SingleElectron/Run2017D-31Mar2018-v1/MINIAOD",
        "isData": True,
        "year": '2017'
    },
    "SingleElectron_Run2017E":{
        "inputDataset": "/SingleElectron/Run2017E-31Mar2018-v1/MINIAOD",
        "isData": True,
        "year": '2017'
    },
    "SingleElectron_Run2017F":{
        "inputDataset": "/SingleElectron/Run2017F-31Mar2018-v1/MINIAOD",
        "isData": True,
        "year": '2017'
    },
    
    
    "EGamma_Run2018A":{
        "inputDataset": "/EGamma/Run2018A-17Sep2018-v2/MINIAOD",
        "isData": True,
        "year": '2018'
    },
    "EGamma_Run2018B":{
        "inputDataset": "/EGamma/Run2018B-17Sep2018-v1/MINIAOD",
        "isData": True,
        "year": '2018'
    },
    "EGamma_Run2018C":{
        "inputDataset": "/EGamma/Run2018C-17Sep2018-v1/MINIAOD",
        "isData": True,
        "year": '2018'
    },
    "EGamma_Run2018D":{
        "inputDataset": "/EGamma/Run2018D-22Jan2019-v2/MINIAOD",
        "isData": True,
        "year": '2018D'
    }
}
    
myJobs = {
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
        "year": '2016'
    },
    "SingleMuon_Run2017B":{
        "inputDataset": "/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD",
        "isData": True,
        "year": '2017'
    },
    "SingleMuon_Run2017C":{
        "inputDataset": "/SingleMuon/Run2017C-31Mar2018-v1/MINIAOD",
        "isData": True,
        "year": '2017'
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
myJobs = {} 
with open("HNL_samples.txt") as f:
    for line in f:
        line = line.rstrip()
        name = line.rsplit('/')[1]+"-2016"
        print(name)
        myJobs[name] = {
            "inputDataset": line,
            "year": '2016',
            "unitsPerJob": 15,
            "isData":False,
            "addLLPInfo": True,
            "addSignalLHE": True
        }
'''


requestName = "NANOX_110320-v2"
userName = "vcepaiti"
#userName = getUsernameFromSiteDB() 
configTmpl = Configuration()

configTmpl.section_('General')
configTmpl.General.transferOutputs = True
configTmpl.General.transferLogs = False

configTmpl.section_('JobType')
configTmpl.JobType.psetName = "LLPReco/NANOProducer/test/produceNANO.py"
configTmpl.JobType.pluginName = 'Analysis'
configTmpl.JobType.outputFiles = ['nano.root']
configTmpl.JobType.allowUndistributedCMSSW = True
configTmpl.JobType.maxJobRuntimeMin= 15*60
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
        print len(jobName+"_"+requestName)
        print config.General.requestName
        config.General.workArea = "crab/"+requestName+"/"+jobName
        config.Data.outLFNDirBase = "/store/user/"+userName+"/LLP/"+requestName+"/"+jobName
        userInputFiles = myJob.get('userInputFiles', None)
        if not userInputFiles:
            config.Data.inputDataset = myJob["inputDataset"]
            if "/USER" not in myJob["inputDataset"]:
                config.Data.inputDBS = myJob.get('inputDBS', 'global')
            else:
                config.Data.inputDBS = "phys03"
            print config.Data.inputDBS
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
            config.Data.unitsPerJob = myJob.get('unitsPerJob', 30)
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
        
