import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("USER")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag (process.GlobalTag, 'auto:run2_mc')

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #"root://cms-xrd-global.cern.ch///store/mc/RunIIFall17MiniAODv2/gluinoGMSB_M2500_ctau1000p0_TuneCP2_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/80000/C65B9A27-9208-E911-8C95-485B39897242.root"
        #"root://cms-xrd-global.cern.ch://store/mc/RunIISummer16MiniAODv2/SMS-T1qqqq_ctau-0p01_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_GridpackScan_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/110000/FE913E88-B389-E711-A5A6-001E67579ED8.root",
        #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mkomm/LLP/LLP_miniaodv180920/SMS-T1qqqq_ctau-10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/SMS-T1qqqq_ctau-10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/LLP_miniaodv180920/180926_191050/0000/T1qqqqLL_MINIAODSIM_82.root'
        'root://cms-xrd-global.cern.ch://store/mc/RunIISummer16MiniAODv2/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/182EF290-92B1-E611-A574-0CC47A7C361E.root',
    )
)

## Output file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.OUT = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    outputCommands = cms.untracked.vstring(['keep *'])
)

from PhysicsTools.PatAlgos.tools.jetTools import *

updateJetCollection(
    process,
    jetSource=cms.InputTag("slimmedJets"),
    jetCorrections=("AK4PFchs", cms.vstring(["L1FastJet", "L2Relative", "L3Absolute"]), "None"),
    btagInfos = ['pfImpactParameterTagInfos', 'pfInclusiveSecondaryVertexFinderTagInfos', 'pfDeepCSVTagInfos']
)

process.updatedPatJets.addBTagInfo = cms.bool(True)
process.updatedPatJets.addDiscriminators = cms.bool(True)
process.updatedPatJets.addJetCorrFactors = cms.bool(True)
process.updatedPatJets.addTagInfos = cms.bool(True)


process.pfXTagInfos = cms.EDProducer("XTagInfoProducer",
    jets = cms.InputTag("updatedPatJets"),
    shallow_tag_infos = cms.InputTag('pfDeepCSVTagInfos'),
    vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
    secondary_vertices = cms.InputTag("slimmedSecondaryVertices")
)

process.pfXTags = cms.EDProducer("XTagProducer",
    graph_path=cms.FileInPath("LLPReco/XTagProducer/data/da.pb"),
    src=cms.InputTag("pfXTagInfos"),
    ctau_values=cms.vdouble(-2., 0., 3.), # provide log(ctau/1mm) to be evaluated: i.e. 10 mum, 1 mm and 1 m here
    ctau_descriptors=cms.vstring("0p01", "1", "1000") # provide log(ctau/1mm) to be evaluated: i.e. 1 mum, 1 mm and 1 m here
)

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
)

process.task = cms.Task(
        process.pfImpactParameterTagInfos,
        process.pfInclusiveSecondaryVertexFinderTagInfos,
        process.pfDeepCSVTagInfos,
        process.patJetCorrFactors,
        process.pfXTagInfos,
        process.pfXTags,
        process.updatedPatJets
)
process.p = cms.Path(process.task)

process.endpath= cms.EndPath(process.OUT)
