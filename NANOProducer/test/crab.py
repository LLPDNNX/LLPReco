import CRABClient
import datetime,sys,os
import copy
import math
import urllib, json
from WMCore.Configuration import Configuration

myJobsTraining = {
    "TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-2016":{
        "inputDataset":"/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
        "year": 2016,
        "unitsPerJob":1,
    },
    "TTJets_TuneCP5_13TeV-madgraphMLM-pythia8-2017":{
        "inputDataset":"/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob":1,
    },
    "TTJets_TuneCP5_13TeV-madgraphMLM-pythia8-2018":{
        "inputDataset":"/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob":1,
    },


    "WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-2016":{
        "inputDataset":"/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob":1,
    },
    "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8-2017":{
        "inputDataset":"/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM",
        "year": 2017,
        "unitsPerJob":1,
    },
    "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8-2018":{
        "inputDataset":"/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
        "year": 2018,
        "unitsPerJob":1,
    },



    "QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob":2,
    },
    "QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob":2,
    },
    "QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8-2016":{
        "inputDataset":"/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "year": 2016,
        "unitsPerJob":2,
    },
    
    
    
    "QCD_Pt_15to30_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt_15to30_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM",
        "year": 2017,
        "unitsPerJob":2,
    },
    "QCD_Pt_30to50_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt_30to50_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM",
        "year": 2017,
        "unitsPerJob":2,
    },
    "QCD_Pt_50to80_TuneCP5_13TeV_pythia8-2017":{
        "inputDataset":"/QCD_Pt_50to80_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM",
        "year": 2017,
        "unitsPerJob":2,
    },
    
    
    
    "QCD_Pt_15to30_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt_15to30_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM",
        "year": 2018,
        "unitsPerJob":2,
    },
    "QCD_Pt_30to50_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt_30to50_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob":2,
    },
    "QCD_Pt_50to80_TuneCP5_13TeV_pythia8-2018":{
        "inputDataset":"/QCD_Pt_50to80_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM",
        "year": 2018,
        "unitsPerJob":2,
    },
    
    
    "GluGluHToTauTau_M125_13TeV_powheg_pythia8-2016": {
        "inputDataset": "/GluGluHToTauTau_M125_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v3/MINIAODSIM",
        "year": 2016,
        "unitsPerJob":1,
    },
    
    "GluGluHToTauTau_M125_13TeV_powheg_pythia8-2017": {
        "inputDataset": "/GluGluHToTauTau_M125_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob":2,
    },
    
    "GluGluHToTauTau_M125_13TeV_powheg_pythia8-2018": {
        "inputDataset": "/GluGluHToTauTau_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
        "year": 2018,
        "unitsPerJob":2,
    },

    
    "LLPGun-2016":{
        "inputDataset":"/LLPGun/mkomm-miniaod16v3_200929-53f8667ba4b240d5eafd36e71bf34742/USER",
        "year": 2016,
        "unitsPerJob": 50,
    },
    "LLPGun-2017":{
        "inputDataset":"/LLPGun/mkomm-miniaod17v2_200929-442a7f6ea2510b243c486adb7160c528/USER",
        "year": 2017,
        "unitsPerJob": 50,
    },
    "LLPGun-2018":{
        "inputDataset":"/LLPGun/mkomm-miniaod18_200929-c21dec93027231dc6f615dfe5c662834/USER",
        "year": 2018,
        "unitsPerJob": 50,
    },
}

myJobsAnalysis = {

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

    "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_ext3-2017":{
        "inputDataset":"/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext3-v1/MINIAODSIM",
        "year": 2017,
        "unitsPerJob": 3
    },

    "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8-2018":{
        "inputDataset":"/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 3
    },

    "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_ext2-2018": {
        "inputDataset":"/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM",
        "year": 2018,
        "unitsPerJob": 3
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
        "year": 2017,
        "unitsPerJob": 2
    },

    "ZGToLLG_01J_5f_PDFWeights_TuneCP5_13TeV-amcatnloFXFX-pythia8-2018":{
        "inputDataset": "/ZGToLLG_01J_5f_PDFWeights_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
        "year": 2018,
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
myJobsData = {
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
myJobsHNL = {} 
with open("HNL_signal2018.txt") as f:
    for line in f:
        line = line.rstrip()
        chunks = line.rsplit('/')
        version = chunks[2]
        if "miniaod16v3" in version:
            year = "2016"
        elif "miniaod17v2" in version:
            year = "2017"
        elif "miniaod18" in version:
            year = "2018"
        name = chunks[1]+"-"+year
        myJobsHNL[name] = {
            "inputDataset": line,
            "year": year,
            "unitsPerJob": 15,
            "isData": False,
            "addLLPInfo": True,
            "addSignalLHE": True
        }
'''
#myJobs = myJobsTraining
#myJobs = myJobsHNL
#myJobs = myJobsData
'''
for k,job in myJobs.items():
    #if str(job['year'])!='2018':
    #    del myJobs[k]
'''


ntau2016 =[
    "/HNL_dirac_ntau0_ctau1p0e-01_massHNL10p0_Vall5p262e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-01_massHNL12p0_Vall3p272e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-01_massHNL14p0_Vall2p193e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-01_massHNL16p0_Vall1p551e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-01_massHNL18p0_Vall1p144e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-01_massHNL20p0_Vall8p709e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-01_massHNL8p0_Vall9p475e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-02_massHNL10p0_Vall1p664e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-02_massHNL12p0_Vall1p035e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-02_massHNL14p0_Vall6p933e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-02_massHNL16p0_Vall4p905e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-02_massHNL18p0_Vall3p617e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-02_massHNL20p0_Vall2p754e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-03_massHNL12p0_Vall3p272e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-03_massHNL14p0_Vall2p193e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-03_massHNL16p0_Vall1p551e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-03_massHNL18p0_Vall1p144e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-03_massHNL20p0_Vall8p709e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-04_massHNL14p0_Vall6p933e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-04_massHNL16p0_Vall4p905e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-04_massHNL18p0_Vall3p617e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e-04_massHNL20p0_Vall2p754e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e00_massHNL10p0_Vall1p664e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e00_massHNL12p0_Vall1p035e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e00_massHNL14p0_Vall6p933e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e00_massHNL16p0_Vall4p905e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e00_massHNL18p0_Vall3p617e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e00_massHNL4p5_Vall1p438e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e00_massHNL6p0_Vall6p496e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e00_massHNL8p0_Vall2p996e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e01_massHNL10p0_Vall5p262e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e01_massHNL12p0_Vall3p272e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e01_massHNL14p0_Vall2p193e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e01_massHNL16p0_Vall1p551e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e01_massHNL18p0_Vall1p144e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e01_massHNL1p0_Vall2p359e-01/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e01_massHNL1p5_Vall8p442e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e01_massHNL2p0_Vall4p066e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e01_massHNL3p0_Vall1p388e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e01_massHNL4p5_Vall4p549e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e01_massHNL6p0_Vall2p054e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e01_massHNL8p0_Vall9p475e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e02_massHNL10p0_Vall1p664e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e02_massHNL1p0_Vall7p460e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e02_massHNL1p5_Vall2p670e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e02_massHNL2p0_Vall1p286e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e02_massHNL3p0_Vall4p390e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e02_massHNL4p5_Vall1p438e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e02_massHNL6p0_Vall6p496e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e02_massHNL8p0_Vall2p996e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e03_massHNL1p0_Vall2p359e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e03_massHNL1p5_Vall8p442e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e03_massHNL2p0_Vall4p066e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e03_massHNL3p0_Vall1p388e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e03_massHNL4p5_Vall4p549e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e03_massHNL6p0_Vall2p054e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e04_massHNL1p0_Vall7p460e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e04_massHNL1p5_Vall2p670e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau0_ctau1p0e04_massHNL2p0_Vall1p286e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-01_massHNL10p0_Vall5p262e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-01_massHNL12p0_Vall3p272e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-01_massHNL14p0_Vall2p193e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-01_massHNL16p0_Vall1p551e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-01_massHNL18p0_Vall1p144e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-01_massHNL20p0_Vall8p709e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-01_massHNL8p0_Vall9p475e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-02_massHNL10p0_Vall1p664e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-02_massHNL12p0_Vall1p035e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-02_massHNL14p0_Vall6p933e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-02_massHNL16p0_Vall4p905e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-02_massHNL18p0_Vall3p617e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-02_massHNL20p0_Vall2p754e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-03_massHNL12p0_Vall3p272e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-03_massHNL14p0_Vall2p193e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-03_massHNL16p0_Vall1p551e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-03_massHNL18p0_Vall1p144e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-03_massHNL20p0_Vall8p709e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-04_massHNL14p0_Vall6p933e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-04_massHNL16p0_Vall4p905e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-04_massHNL18p0_Vall3p617e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e-04_massHNL20p0_Vall2p754e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e00_massHNL10p0_Vall1p664e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e00_massHNL12p0_Vall1p035e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e00_massHNL14p0_Vall6p933e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e00_massHNL16p0_Vall4p905e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e00_massHNL18p0_Vall3p617e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e00_massHNL4p5_Vall1p438e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e00_massHNL6p0_Vall6p496e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e00_massHNL8p0_Vall2p996e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e01_massHNL10p0_Vall5p262e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e01_massHNL12p0_Vall3p272e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e01_massHNL14p0_Vall2p193e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e01_massHNL16p0_Vall1p551e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e01_massHNL18p0_Vall1p144e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e01_massHNL1p0_Vall2p359e-01/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e01_massHNL1p5_Vall8p442e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e01_massHNL2p0_Vall4p066e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e01_massHNL3p0_Vall1p388e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e01_massHNL4p5_Vall4p549e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e01_massHNL6p0_Vall2p054e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e01_massHNL8p0_Vall9p475e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e02_massHNL10p0_Vall1p664e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e02_massHNL1p0_Vall7p460e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e02_massHNL1p5_Vall2p670e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e02_massHNL2p0_Vall1p286e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e02_massHNL3p0_Vall4p390e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e02_massHNL4p5_Vall1p438e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e02_massHNL6p0_Vall6p496e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e02_massHNL8p0_Vall2p996e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e03_massHNL1p0_Vall2p359e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e03_massHNL1p5_Vall8p442e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e03_massHNL2p0_Vall4p066e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e03_massHNL3p0_Vall1p388e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e03_massHNL4p5_Vall4p549e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e03_massHNL6p0_Vall2p054e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e04_massHNL1p0_Vall7p460e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e04_massHNL1p5_Vall2p670e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau1_ctau1p0e04_massHNL2p0_Vall1p286e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-01_massHNL10p0_Vall5p262e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-01_massHNL12p0_Vall3p272e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-01_massHNL14p0_Vall2p193e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-01_massHNL16p0_Vall1p551e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-01_massHNL18p0_Vall1p144e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-01_massHNL20p0_Vall8p709e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-01_massHNL8p0_Vall9p475e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-02_massHNL10p0_Vall1p664e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-02_massHNL12p0_Vall1p035e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-02_massHNL14p0_Vall6p933e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-02_massHNL16p0_Vall4p905e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-02_massHNL18p0_Vall3p617e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-02_massHNL20p0_Vall2p754e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-03_massHNL12p0_Vall3p272e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-03_massHNL14p0_Vall2p193e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-03_massHNL16p0_Vall1p551e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-03_massHNL18p0_Vall1p144e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-03_massHNL20p0_Vall8p709e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-04_massHNL14p0_Vall6p933e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-04_massHNL16p0_Vall4p905e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-04_massHNL18p0_Vall3p617e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e-04_massHNL20p0_Vall2p754e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e00_massHNL10p0_Vall1p664e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e00_massHNL12p0_Vall1p035e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e00_massHNL14p0_Vall6p933e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e00_massHNL16p0_Vall4p905e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e00_massHNL18p0_Vall3p617e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e00_massHNL4p5_Vall1p438e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e00_massHNL6p0_Vall6p496e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e00_massHNL8p0_Vall2p996e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e01_massHNL10p0_Vall5p262e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e01_massHNL12p0_Vall3p272e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e01_massHNL14p0_Vall2p193e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e01_massHNL16p0_Vall1p551e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e01_massHNL18p0_Vall1p144e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e01_massHNL1p0_Vall2p359e-01/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e01_massHNL1p5_Vall8p442e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e01_massHNL3p0_Vall1p388e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e01_massHNL4p5_Vall4p549e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e01_massHNL6p0_Vall2p054e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e01_massHNL8p0_Vall9p475e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e02_massHNL10p0_Vall1p664e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e02_massHNL1p0_Vall7p460e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e02_massHNL2p0_Vall1p286e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e02_massHNL3p0_Vall4p390e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e02_massHNL4p5_Vall1p438e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e02_massHNL6p0_Vall6p496e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e02_massHNL8p0_Vall2p996e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e03_massHNL3p0_Vall1p388e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e03_massHNL4p5_Vall4p549e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_dirac_ntau2_ctau1p0e03_massHNL6p0_Vall2p054e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-01_massHNL10p0_Vall3p721e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-01_massHNL12p0_Vall2p314e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-01_massHNL14p0_Vall1p551e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-01_massHNL16p0_Vall1p097e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-01_massHNL18p0_Vall8p092e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-01_massHNL20p0_Vall6p165e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-01_massHNL8p0_Vall6p702e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-02_massHNL10p0_Vall1p177e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-02_massHNL12p0_Vall7p319e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-02_massHNL14p0_Vall4p904e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-02_massHNL16p0_Vall3p470e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-02_massHNL18p0_Vall2p559e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-02_massHNL20p0_Vall1p950e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-03_massHNL12p0_Vall2p314e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-03_massHNL14p0_Vall1p551e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-03_massHNL16p0_Vall1p097e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-03_massHNL18p0_Vall8p092e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-03_massHNL20p0_Vall6p165e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-04_massHNL14p0_Vall4p904e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-04_massHNL16p0_Vall3p470e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-04_massHNL18p0_Vall2p559e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e-04_massHNL20p0_Vall1p950e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e00_massHNL10p0_Vall1p177e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e00_massHNL12p0_Vall7p319e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e00_massHNL14p0_Vall4p904e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e00_massHNL16p0_Vall3p470e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e00_massHNL18p0_Vall2p559e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e00_massHNL4p5_Vall1p016e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e00_massHNL6p0_Vall4p597e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e00_massHNL8p0_Vall2p119e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e01_massHNL10p0_Vall3p721e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e01_massHNL12p0_Vall2p314e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e01_massHNL14p0_Vall1p551e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e01_massHNL16p0_Vall1p097e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e01_massHNL18p0_Vall8p092e-05/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e01_massHNL1p0_Vall1p668e-01/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e01_massHNL1p5_Vall5p965e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e01_massHNL2p0_Vall2p871e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e01_massHNL3p0_Vall9p825e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e01_massHNL4p5_Vall3p213e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e01_massHNL6p0_Vall1p454e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e01_massHNL8p0_Vall6p702e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e02_massHNL10p0_Vall1p177e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e02_massHNL1p0_Vall5p274e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e02_massHNL1p5_Vall1p886e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e02_massHNL2p0_Vall9p078e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e02_massHNL3p0_Vall3p107e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e02_massHNL4p5_Vall1p016e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e02_massHNL6p0_Vall4p597e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e02_massHNL8p0_Vall2p119e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e03_massHNL1p0_Vall1p668e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e03_massHNL1p5_Vall5p965e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e03_massHNL2p0_Vall2p871e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e03_massHNL3p0_Vall9p825e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e03_massHNL4p5_Vall3p213e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e03_massHNL6p0_Vall1p454e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e04_massHNL1p0_Vall5p274e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e04_massHNL1p5_Vall1p886e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau0_ctau1p0e04_massHNL2p0_Vall9p078e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-01_massHNL10p0_Vall3p721e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-01_massHNL12p0_Vall2p314e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-01_massHNL14p0_Vall1p551e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-01_massHNL16p0_Vall1p097e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-01_massHNL18p0_Vall8p092e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-01_massHNL20p0_Vall6p165e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-01_massHNL8p0_Vall6p702e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-02_massHNL10p0_Vall1p177e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-02_massHNL12p0_Vall7p319e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-02_massHNL14p0_Vall4p904e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-02_massHNL16p0_Vall3p470e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-02_massHNL18p0_Vall2p559e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-02_massHNL20p0_Vall1p950e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-03_massHNL12p0_Vall2p314e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-03_massHNL14p0_Vall1p551e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-03_massHNL16p0_Vall1p097e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-03_massHNL18p0_Vall8p092e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-03_massHNL20p0_Vall6p165e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-04_massHNL14p0_Vall4p904e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-04_massHNL16p0_Vall3p470e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-04_massHNL18p0_Vall2p559e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e-04_massHNL20p0_Vall1p950e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e00_massHNL10p0_Vall1p177e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e00_massHNL12p0_Vall7p319e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e00_massHNL14p0_Vall4p904e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e00_massHNL16p0_Vall3p470e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e00_massHNL18p0_Vall2p559e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e00_massHNL4p5_Vall1p016e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e00_massHNL6p0_Vall4p597e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e00_massHNL8p0_Vall2p119e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e01_massHNL10p0_Vall3p721e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e01_massHNL12p0_Vall2p314e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e01_massHNL14p0_Vall1p551e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e01_massHNL16p0_Vall1p097e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e01_massHNL18p0_Vall8p092e-05/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e01_massHNL1p0_Vall1p668e-01/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e01_massHNL1p5_Vall5p965e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e01_massHNL2p0_Vall2p871e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e01_massHNL3p0_Vall9p825e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e01_massHNL4p5_Vall3p213e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e01_massHNL6p0_Vall1p454e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e01_massHNL8p0_Vall6p702e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e02_massHNL10p0_Vall1p177e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e02_massHNL1p0_Vall5p274e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e02_massHNL1p5_Vall1p886e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e02_massHNL2p0_Vall9p078e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e02_massHNL3p0_Vall3p107e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e02_massHNL4p5_Vall1p016e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e02_massHNL6p0_Vall4p597e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e02_massHNL8p0_Vall2p119e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e03_massHNL1p0_Vall1p668e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e03_massHNL1p5_Vall5p965e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e03_massHNL2p0_Vall2p871e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e03_massHNL3p0_Vall9p825e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e03_massHNL4p5_Vall3p213e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e03_massHNL6p0_Vall1p454e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e04_massHNL1p0_Vall5p274e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e04_massHNL1p5_Vall1p886e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau1_ctau1p0e04_massHNL2p0_Vall9p078e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-01_massHNL10p0_Vall3p721e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-01_massHNL12p0_Vall2p314e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-01_massHNL14p0_Vall1p551e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-01_massHNL16p0_Vall1p097e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-01_massHNL18p0_Vall8p092e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-01_massHNL20p0_Vall6p165e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-01_massHNL8p0_Vall6p702e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-02_massHNL10p0_Vall1p177e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-02_massHNL12p0_Vall7p319e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-02_massHNL14p0_Vall4p904e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-02_massHNL16p0_Vall3p470e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-02_massHNL18p0_Vall2p559e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-02_massHNL20p0_Vall1p950e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-03_massHNL12p0_Vall2p314e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-03_massHNL14p0_Vall1p551e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-03_massHNL16p0_Vall1p097e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-03_massHNL18p0_Vall8p092e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-03_massHNL20p0_Vall6p165e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-04_massHNL14p0_Vall4p904e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-04_massHNL16p0_Vall3p470e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-04_massHNL18p0_Vall2p559e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e-04_massHNL20p0_Vall1p950e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e00_massHNL10p0_Vall1p177e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e00_massHNL12p0_Vall7p319e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e00_massHNL12p0_Vall7p319e-04/mkomm-miniaod16v3_221221-dbc7933c2a7a8185653b951ff340c9d2/USER",
    "/HNL_majorana_ntau2_ctau1p0e00_massHNL14p0_Vall4p904e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e00_massHNL16p0_Vall3p470e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e00_massHNL18p0_Vall2p559e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e00_massHNL4p5_Vall1p016e-02/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e00_massHNL6p0_Vall4p597e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e00_massHNL8p0_Vall2p119e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e01_massHNL10p0_Vall3p721e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e01_massHNL12p0_Vall2p314e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e01_massHNL14p0_Vall1p551e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e01_massHNL16p0_Vall1p097e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e01_massHNL18p0_Vall8p092e-05/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e01_massHNL3p0_Vall9p825e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e01_massHNL4p5_Vall3p213e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e01_massHNL6p0_Vall1p454e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e01_massHNL8p0_Vall6p702e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e02_massHNL10p0_Vall1p177e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e02_massHNL3p0_Vall3p107e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e02_massHNL4p5_Vall1p016e-03/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e02_massHNL6p0_Vall4p597e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e02_massHNL8p0_Vall2p119e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e03_massHNL3p0_Vall9p825e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e03_massHNL4p5_Vall3p213e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER",
    "/HNL_majorana_ntau2_ctau1p0e03_massHNL6p0_Vall1p454e-04/mkomm-miniaod16v3_221221-53f8667ba4b240d5eafd36e71bf34742/USER"
]

myJobs = {}
for ntauDataset in ntau2016:
    sampleName = ntauDataset.split('/')[1]+"-2016"
    #print sampleName
    myJobs[sampleName] = {
        "inputDataset":ntauDataset,
        "year": 2016,
        "unitsPerJob": 50,
        "addSignalLHE": True
    }
    
requestName = "NANOX_221221"
userName = "mkomm"
configTmpl = Configuration()

configTmpl.section_('General')
configTmpl.General.transferOutputs = True
configTmpl.General.transferLogs = False

configTmpl.section_('JobType')
configTmpl.JobType.psetName = "LLPReco/NANOProducer/test/produceNANO.py"
configTmpl.JobType.pluginName = 'Analysis'
configTmpl.JobType.outputFiles = ['nano.root']
configTmpl.JobType.allowUndistributedCMSSW = True
configTmpl.JobType.pyCfgParams = []
configTmpl.JobType.inputFiles = []
configTmpl.JobType.maxMemoryMB = 2500
configTmpl.JobType.maxJobRuntimeMin= 10*60
configTmpl.section_('Data')
configTmpl.Data.useParent = False
configTmpl.section_('Site')
configTmpl.Site.storageSite = 'T2_UK_London_IC'

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
        if i>400:
            break
        print(jobName)
        isData = False
        myJob = myJobs[jobName]
        i=i+1
        config = copy.deepcopy(configTmpl)
        config.General.requestName = jobName
        config.General.workArea = "crab/"+requestName+"/"+jobName
        config.Data.outLFNDirBase = "/store/user/"+userName+"/HNL/"+requestName+"/"+jobName
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
        addSignalLHE = str(myJob.get('addSignalLHE', False))
        
        if year not in ['2016','2017','2018','2018D']:
            print "ERROR: Year invalid: ", year
            continue
        print "year:", year
        
        config.JobType.pyCfgParams.append("year="+str(year))
        config.JobType.pyCfgParams.append("addSignalLHE={}".format(addSignalLHE))

        if isData:
            config.JobType.pyCfgParams.append("isData=True")
            config.Data.splitting = 'LumiBased'
            config.Data.unitsPerJob = myJob.get('unitsPerJob', 30)
            config.JobType.maxJobRuntimeMin = myJob.get('maxJobRuntimeMin', 16*60)
            if year == '2016':
                config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification//Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
            if year == '2017':
                config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
            if year == '2018' or year== '2018D':
                config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt' 

            # Recovery task!
            #config.Data.lumiMask = '{}.json'.format(jobName) 

        else:
            config.JobType.pyCfgParams.append("isData=False")
            if myJob.has_key('unitsPerJob'):
                config.Data.splitting = 'FileBased'
                config.JobType.maxJobRuntimeMin= 10*60
                config.Data.unitsPerJob = myJob.get('unitsPerJob', 1)
            else:
                config.Data.splitting = 'Automatic'
                config.Data.unitsPerJob = 10*60

        if "params" in myJob:
            params = myJob["params"]

            for param in params:
                config.JobType.pyCfgParams.append(str(param))
                
        if "whitelist" in myJob:
            config.Site.whitelist = myJob['whitelist']
            
        if "blacklist" in myJob:
            config.Site.blacklist = myJob['blacklist']

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
        
        #break
        
        print
        print
        
