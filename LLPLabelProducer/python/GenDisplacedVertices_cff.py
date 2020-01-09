import FWCore.ParameterSet.Config as cms


#merge GenParticles for easy matching of GenJets to vertices
genParticlesMerged = cms.EDProducer("MergedGenParticleProducer",
    inputPruned = cms.InputTag("prunedGenParticles"),
    inputPacked = cms.InputTag("packedGenParticles")
)


#select only stable particles (status=1) for reclustering GenJets for the matching; exclude neutrinos & SUSY particles
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
genParticlesForGenJets = cms.EDFilter(
    "CandPtrSelector", 
    src = cms.InputTag("genParticlesMerged"),
    cut = cms.string("status==1 && abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16 && abs(pdgId) < 1000000"),
    jetPtMin = cms.double(3.0)
)


#recluster GenJets
genJetsReclustered = ak4GenJets.clone(
    src = cms.InputTag('genParticlesForGenJets')
)

#produce DisplacedGenVertices and match to GenJets
displacedGenVertices = cms.EDProducer(
    "DisplacedGenVertexProducer",
    srcGenParticles = cms.InputTag("genParticlesMerged"),
    srcGenJets = cms.InputTag("genJetsReclustered")
)

displacedGenVertexSequence = cms.Sequence(
    genParticlesMerged     +
    genParticlesForGenJets +
    genJetsReclustered     +
    displacedGenVertices
)
