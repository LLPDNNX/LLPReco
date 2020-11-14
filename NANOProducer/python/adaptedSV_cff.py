import FWCore.ParameterSet.Config as cms

# https://gist.github.com/daiktas/d3e0ed47b82e45502aaf787ae12934fd

adaptedInclusiveCandidateVertexFinder = cms.EDProducer("InclusiveCandidateVertexFinder",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterizer = cms.PSet(
        clusterMaxDistance = cms.double(0.05),
        clusterMaxSignificance = cms.double(4.5),
        clusterMinAngleCosine = cms.double(0.5),
        distanceRatio = cms.double(20),
        seedMax3DIPSignificance = cms.double(9999.0),
        seedMax3DIPValue = cms.double(9999.0),
        seedMin3DIPSignificance = cms.double(1.2),
        seedMin3DIPValue = cms.double(0.005)
    ),
    fitterRatio = cms.double(0.25),
    fitterSigmacut = cms.double(3),
    fitterTini = cms.double(256),
    maxNTracks = cms.uint32(30),
    maximumLongitudinalImpactParameter = cms.double(1), # Default is 0.3->20 (Ghent) -> 1 may be safer
    minHits = cms.uint32(6), # 8 -> 6
    minPt = cms.double(0.8),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    tracks = cms.InputTag("packedPFCandidates"),
    useDirectVertexFitter = cms.bool(True),
    useVertexReco = cms.bool(True),
    vertexMinAngleCosine = cms.double(0.95),
    vertexMinDLen2DSig = cms.double(2.5),
    vertexMinDLenSig = cms.double(0.5),
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        primcut = cms.double(1.0),
        seccut = cms.double(3),
        smoothing = cms.bool(True)
    )
)

adaptedCandidateVertexMerger = cms.EDProducer("CandidateVertexMerger",
    maxFraction = cms.double(0.7),
    minSignificance = cms.double(2),
    secondaryVertices = cms.InputTag("adaptedInclusiveCandidateVertexFinder")
)

adaptedCandidateVertexArbitrator = cms.EDProducer("CandidateVertexArbitrator",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    dLenFraction = cms.double(0.333),
    dRCut = cms.double(0.4), # 0.4 -> 1 # not sure good idea...
    distCut = cms.double(0.04), # 0.04 -> 0.1 # check if safe
    fitterRatio = cms.double(0.25),
    fitterSigmacut = cms.double(3),
    fitterTini = cms.double(256),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    secondaryVertices = cms.InputTag("adaptedCandidateVertexMerger"),
    sigCut = cms.double(5),
    trackMinLayers = cms.int32(4),
    trackMinPixels = cms.int32(0), # 1->0
    trackMinPt = cms.double(0.8),
    tracks = cms.InputTag("packedPFCandidates")
)

adaptedSlimmedSecondaryVertices = cms.EDProducer("CandidateVertexMerger",
    maxFraction = cms.double(0.2),
    minSignificance = cms.double(10.0),
    secondaryVertices = cms.InputTag("adaptedCandidateVertexArbitrator")
)

adaptedVertexing = cms.Sequence(adaptedInclusiveCandidateVertexFinder +
    adaptedCandidateVertexMerger +
    adaptedCandidateVertexArbitrator +
    adaptedSlimmedSecondaryVertices
)