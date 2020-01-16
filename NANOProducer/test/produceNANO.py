import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
#from PhysicsTools.NanoAOD.common_cff import *
from Configuration.StandardSequences.Eras import eras
from RecoBTag.Configuration.RecoBTag_cff import *
from Configuration.AlCa.GlobalTag import GlobalTag


options = VarParsing ('analysis')

options.register(
    'isData',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "is data"
)

options.register(
    'addLLPInfo',
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "add llp info"
)

options.register(
    'LLPtype',
    "HNL",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "add llp type"
)
options.register(
    'year',
    2017,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "add year file"
)
options.parseArguments() 

if options.year == 2016:
    print "2016"
    process = cms.Process('NANO',eras.Run2_2016,eras.run2_nanoAOD_94X2016)
elif options.year == 2017:
    print "2017"
    process = cms.Process('NANO',eras.Run2_2017,eras.run2_nanoAOD_94XMiniAODv2)
elif options.year == 2018 or options.year == 2019:
    print "2018D"
    process = cms.Process('NANO',eras.Run2_2018,eras.run2_nanoAOD_102Xv1)


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

if options.isData:
    process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
    dataTier = cms.untracked.string('NANOAOD')
else:
    process.load('SimGeneral.MixingModule.mixNoPU_cfi')
    process.load('Configuration.StandardSequences.MagneticField_cff')
    dataTier = cms.untracked.string('NANOAODSIM')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

## Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'root://maite.iihe.ac.be///store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-1_V-0.13416407865_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_91.root',
        #'root://maite.iihe.ac.be///store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-8_V-0.000415932686862_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_140.root',
        'root://maite.iihe.ac.be///store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-20_V-1.22474487139e-05_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_56.root'
    )
)
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('test102X nevts:1000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

process.plain=cms.Path()

def addModule(m):
    process.plain+=m

'''
if not options.isData:
    process.TFileService = cms.Service("TFileService",
        fileName = cms.string("info.root")
    )

    process.eventAndPuInfo = cms.EDAnalyzer("EventInfoCollector",
        GenEventInfo=cms.InputTag("generator"),
        PileupSummaryInfo=cms.InputTag("slimmedAddPileupInfo")
    )
    process.plain+=process.eventAndPuInfo

'''
#Output definition
process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    saveProvenance = cms.untracked.bool(True),
    fakeNameForCrab = cms.untracked.bool(True),
    dataset = cms.untracked.PSet(
        dataTier = dataTier,
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('nano.root'),
    #outputCommands = process.NANOAODSIMEventContent.outputCommands+cms.untracked.vstring(
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep nanoaodFlatTable_*Table_*_*',
        'keep edmTriggerResults_*_*_*',
        'keep nanoaodMergeableCounterTable_*Table_*_*',
        'keep nanoaodUniqueString_nanoMetadata_*_*',
        
        'drop *_caloMetTable_*_*',
        'drop *_fatJetTable_*_*',
        'drop *_genJetAK8FlavourTable_*_*',
        'drop *_genJetAK8Table_*_*',
        'drop *_genVisTauTable_*_*',
        'drop *_subJetTable_*_*',
        'drop *_tkMetTable_*_*',
        'drop *_puppiMetTable_*_*',
        'drop *_ttbarCategoryTable_*_*',
        
        'drop *_photonTable_*_*',
        'drop *_photonMCTable_*_*',
        
        'drop *_tauTable_*_*',
        'drop *_tauMCTable_*_*' ,
        
        'drop *_saJetTable_*_*',
        'drop *_FatJetTable_*_*',
        'drop *_saTable_*_*',
        
        'drop *_simpleCleanerTable_photons_*',
        'drop *_simpleCleanerTable_taus_*',
        
        'drop *_rivetMetTable_*_*',
        'drop *_rivetLeptonTable_*_*',
        'drop *_rivetProducerHTXS_*_*'
        
    )
)

process.MINIAODoutput = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('output.root'),
    outputCommands = process.NANOAODSIMoutput.outputCommands,
    dropMetaData = cms.untracked.string('ALL'),
)


## Output file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.OUT = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    outputCommands = cms.untracked.vstring(['keep *'])
)

if options.isData:
    if options.year == 2016:
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v11', '')
    if options.year == 2017:
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v11', '')
    if options.year == 2018:
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v12', '')
    if options.year == 2019:
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Prompt_v15', '')

    jetCorrectionsAK4PFchs = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'], 'None')
    jetCorrectionsAK4PF = ('AK4PF', ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'], 'None')
else:
    if options.year == 2016:
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_mcRun2_asymptotic_v7', '')
    if options.year == 2017:
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_mc2017_realistic_v7', '')
    if options.year == 2018 or options.year == 2019:
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v20', '')

    jetCorrectionsAK4PFchs = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
    jetCorrectionsAK4PF = ('AK4PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
    

from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.NanoAOD.common_cff import *

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
    muonSrc  = cms.InputTag("slimmedMuons"),
    electronSrc = cms.InputTag("slimmedElectrons"),
    shallow_tag_infos = cms.InputTag('pfDeepCSVTagInfos'),
    vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
    secondary_vertices = cms.InputTag("slimmedSecondaryVertices")
)

process.pfXTags = cms.EDProducer("XTagProducer",
    graph_path=cms.FileInPath("LLPReco/XTagProducer/data/da.pb"),
    src=cms.InputTag("pfXTagInfos"),
    ctau_values=cms.vdouble(-2., 0., 3.), # provide log(ctau/1mm) to be evaluated: i.e. 10 mum, 1 mm and 1 m here
    ctau_descriptors=cms.vstring("0p01", "1", "1000") # provide log(ctau/1mm) to be evaluated: i.e. 10 mum, 1 mm and 1 m here
)

process.nanoTable = cms.EDProducer("NANOProducer",
    srcTags = cms.InputTag("pfXTagInfos"),
    srcLabels = cms.InputTag("llpLabels")
)

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
)

#remove unneeded modules
for moduleName in [
    "genJetAK8Table",
    "genJetAK8FlavourAssociation",
    "genJetAK8FlavourTable",
   
    "particleLevel",
    "rivetLeptonTable",
    "HTXSCategoryTable",
    "rivetMetTable",
    "rivetLeptonTable",
    "rivetProducerHTXS",
    "tautagger",
    "genSubJetAK8Table",
]:
    if hasattr(process,moduleName):
        print "removing module: ",moduleName
        process.nanoSequence.remove(getattr(process,moduleName))
        process.nanoSequenceMC.remove(getattr(process,moduleName))

process.load('LLPReco.LLPLabelProducer.GenDisplacedVertices_cff')

process.llpGenDecayInfo = cms.EDProducer(
    "LLPGenDecayInfoProducer",
    src = cms.InputTag("genParticlesMerged"),
    decays = cms.PSet(
        #hnl -> qql
        hnl = cms.PSet(
            llpId = cms.int32(9990012),
            daughterIds = cms.vint32([1,2,3,4,5,11,13])
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
            daughterIds = cms.vint32([5,11,13])
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
    srcFlavourInfo = cms.InputTag("llpFlavour")
)




if options.isData:
    process.nanoAOD_step = cms.Path(
        process.patJetCorrFactors+
        process.updatedPatJets+
        process.pfXTagInfos+
        process.nanoTable+
        process.nanoSequence
    )
else:
    process.nanoAOD_step = cms.Path(
        process.patJetCorrFactors+
        process.updatedPatJets+
        process.pfXTagInfos+
        process.displacedGenVertexSequence+
        process.llpGenDecayInfo+
        process.llpFlavour+
        process.llpLabels+
        process.nanoTable+
        process.nanoSequenceMC
    )
    
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

process.schedule = cms.Schedule(process.plain,process.nanoAOD_step,process.endjob_step,process.NANOAODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeData,nanoAOD_customizeMC

#call to customisation function nanoAOD_customizeMC imported from PhysicsTools.NanoAOD.nano_cff
if options.isData:
    process = nanoAOD_customizeData(process)
else:
    process = nanoAOD_customizeMC(process)

# End of customisation functions

# Customisation from command line

process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
#print process.dumpPython()
# End adding early deletion

#process.endpath= cms.EndPath(process.OUT)
