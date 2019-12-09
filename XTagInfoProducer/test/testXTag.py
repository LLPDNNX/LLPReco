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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

## Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-15_V-6.76017751246e-05_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_9.root',
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
    graph_path=cms.FileInPath("LLPReco/TensorFlow/data/da.pb"),
    src=cms.InputTag("pfXTagInfos")
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
