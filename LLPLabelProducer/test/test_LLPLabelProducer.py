import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.inputFiles = 'root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-10_V-0.00107238052948_e_Dirac_massiveAndCKM_LO/heavyNeutrino_10.root'
options.parseArguments()

process = cms.Process("TEST")

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        "https://github.com/LLPDNNX/test-files/raw/master/miniaod/Moriond17_aug2018_miniAODv3_HNL.root"
        #"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/SMS-T1qqqq_ctau-1_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_GridpackScan_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/10000/023CA233-4289-E711-9B01-002590FD5A48.root"
        #"/store/user/lbenato/GluGluH_HToSSTobbbb_MH-125_MS-50_ctauS-5_Summer16_MINIAODSIM_calojets/GluGluH_HToSSTobbbb_MH-125_MS-50_ctauS-5_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets/181203_133916/0000/miniaod_6.root"
        #"/store/user/lbenato/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-1_Summer16_MINIAODSIM_calojets/GluGluH_HToSSTobbbb_MH-125_MS-15_ctauS-1_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/RunIISummer16-PU_premix-Moriond17_80X_mcRun2_2016_MINIAODSIM_calojets/181203_134848/0000/miniaod_1.root"
        #"root://maite.iihe.ac.be///store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-1_V-0.13416407865_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_79.root"
        #'root://maite.iihe.ac.be///store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-8_V-0.000415932686862_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_78.root',
        #'/store/mc/RunIISummer16MiniAODv2/DisplacedSUSY_StopToBL_M-1200_CTau-1000_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/72CEA543-D1D5-E611-9956-02163E019CFB.root'
        #'/store/mc/RunIISummer16MiniAODv2/gluinoGMSB_M2500_ctau1p0_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/5A9F1E8C-B6CB-E611-9762-B083FED42C03.root'
    )   
)

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True)
)

process.load('LLPReco.LLPLabelProducer.GenDisplacedVertices_cff')

process.llpGenDecayInfo = cms.EDProducer(
    "LLPGenDecayInfoProducer",
    src = cms.InputTag("genParticlesMerged"),
    decays = cms.PSet(
        #hnl -> qql
        hnl = cms.PSet(
            llpId = cms.int32(9990012),
            daughterIds = cms.vint32([1,2,3,4,5,11,13,15])
        ),
        #gluino -> qq chi0
        split = cms.PSet(
            llpId = cms.int32(1000021),
            daughterIds = cms.vint32([1,2,3,4,5])
        ),
        #gluino -> g gravitino
        gmsb = cms.PSet(
            llpId = cms.int32(1000021),
            daughterIds = cms.vint32([21])
        ),
        #stop -> bl
        rpv = cms.PSet(
            llpId = cms.int32(1000006),
            daughterIds = cms.vint32([5,11,13,15])
        ),
        #H->SS->bbbb
        hss = cms.PSet(
            llpId = cms.int32(9000006),
            daughterIds = cms.vint32([5])
        ),
    )
)

process.llpFlavour = cms.EDProducer(
    "LLPGhostFlavourProducer",
    srcJets = cms.InputTag("slimmedJets"),
    srcDecayInfo = cms.InputTag("llpGenDecayInfo"),
    jetAlgorithm = cms.string("AntiKt"),
    rParam = cms.double(0.4),
    ghostRescaling = cms.double(1e-18),
    relPtTolerance = cms.double(1e-3)
)

process.llpLabels = cms.EDProducer(
    "LLPLabelProducer",
    srcVertices = cms.InputTag("displacedGenVertices"),
    srcJets = cms.InputTag("slimmedJets"),
    srcFlavourInfo = cms.InputTag("llpFlavour"),
    tauPtThreshold = cms.double(1.),
    quarkPtThreshold = cms.double(1.),
    bPtThreshold = cms.double(1.),
    muonPtThreshold = cms.double(1.),
    electronPtThreshold = cms.double(1.),
)


process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",                     
    compressionAlgorithm = cms.untracked.string('LZMA'),    
    compressionLevel = cms.untracked.int32(4),              
    dataset = cms.untracked.PSet(   
        dataTier = cms.untracked.string(''),                
        filterName = cms.untracked.string('')               
    ),      
    dropMetaData = cms.untracked.string('ALL'),             
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),                   
    fastCloning = cms.untracked.bool(False),                
    fileName = cms.untracked.string('test.root'),  
    outputCommands = cms.untracked.vstring(
       'drop *',
       'keep *_*_*_TEST',
    ),
    overrideInputFileSplitLevels = cms.untracked.bool(True) 
)           

process.path = cms.Path(
    process.displacedGenVertexSequence
    *process.llpGenDecayInfo
    *process.llpFlavour
    *process.llpLabels    
) 
process.endpath = cms.EndPath(process.MINIAODSIMoutput) 

           
