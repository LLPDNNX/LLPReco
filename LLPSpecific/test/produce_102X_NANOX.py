import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
#from PhysicsTools.NanoAOD.common_cff import *
from Configuration.StandardSequences.Eras import eras
from RecoBTag.Configuration.RecoBTag_cff import *

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
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

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

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring()
)

process.source.fileNames = [
        'root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-15_V-6.76017751246e-05_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_1.root',
        'root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-15_V-6.76017751246e-05_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_2.root',
        'root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-15_V-6.76017751246e-05_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_3.root',
        'root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-15_V-6.76017751246e-05_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_4.root',
        'root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-15_V-6.76017751246e-05_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_5.root',
        'root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-15_V-6.76017751246e-05_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_6.root',
        'root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-15_V-6.76017751246e-05_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_7.root',
        'root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-15_V-6.76017751246e-05_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_8.root',
        'root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-15_V-6.76017751246e-05_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_9.root',
        'root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-15_V-6.76017751246e-05_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_10.root'
        #'root://cms-xrd-global.cern.ch///store/mc/RunIIFall17MiniAODv2/gluinoGMSB_M2500_ctau1000p0_TuneCP2_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/80000/C65B9A27-9208-E911-8C95-485B39897242.root'
        #'root://cms-xrd-global.cern.ch///store/mc/RunIIFall17MiniAODv2/gluinoGMSB_M2500_ctau1000p0_TuneCP2_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/10D10064-6302-E911-8F33-6C3BE5B58198.root'
        #'root://cms-xrd-global.cern.ch///store/data/Run2017E/SingleMuon/MINIAOD/31Mar2018-v1/00000/947DBA65-773A-E811-A9F8-0023AEEEB703.root'
        ]

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('test102X nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

process.plain=cms.Path()

def addModule(m):
    process.plain+=m

if not options.isData:
    process.TFileService = cms.Service("TFileService", 
        fileName = cms.string("info.root")
    )

    process.eventAndPuInfo = cms.EDAnalyzer("EventInfoCollector",
        GenEventInfo=cms.InputTag("generator"),
        PileupSummaryInfo=cms.InputTag("slimmedAddPileupInfo")
    )
    process.plain+=process.eventAndPuInfo

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

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag

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
    
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    
bTagInfos = [
    'pfImpactParameterTagInfos',
    'pfInclusiveSecondaryVertexFinderTagInfos',
    'pfDeepCSVTagInfos',
    #'pfDeepFlavourTagInfos',
]
bTagDiscriminators = [
    'pfCombinedSecondaryVertexV2BJetTags',
    'pfDeepCSVJetTags:probb',
    'pfDeepCSVJetTags:probc',
    'pfDeepCSVJetTags:probbb',
    #'pfDeepFlavourJetTags:probb',
    #'pfDeepFlavourJetTags:probbb',
    #'pfDeepFlavourJetTags:problepb',
]


updateJetCollection(
        process,
        labelName = "XTag",
        jetSource = cms.InputTag('slimmedJets'),#'ak4Jets'
        jetCorrections = jetCorrectionsAK4PFchs,
        pfCandidates = cms.InputTag('packedPFCandidates'),
        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
)

#there seems to be a bug: addTagInfos is set to false despite len(btagInfos)>0
#to be sure specify included information explicitly here
process.updatedPatJetsTransientCorrectedXTag.addBTagInfo = cms.bool(True)
process.updatedPatJetsTransientCorrectedXTag.addDiscriminators = cms.bool(True)
process.updatedPatJetsTransientCorrectedXTag.addJetCorrFactors = cms.bool(True)
process.updatedPatJetsTransientCorrectedXTag.addTagInfos = cms.bool(True)


process.updateJetXTagSequence = cms.Sequence(
    process.patJetCorrFactorsXTag
    +process.updatedPatJetsXTag
    +process.pfImpactParameterTagInfosXTag
    +process.pfInclusiveSecondaryVertexFinderTagInfosXTag
    +process.pfDeepCSVTagInfosXTag
    +process.pfDeepCSVJetTagsXTag
    +process.pfCombinedSecondaryVertexV2BJetTagsXTag
    +process.patJetCorrFactorsTransientCorrectedXTag
    #+process.pfDeepFlavourTagInfosXTag
    #+process.pfDeepFlavourJetTagsXTag
    +process.updatedPatJetsTransientCorrectedXTag
)

addModule(process.updateJetXTagSequence)


process.selectJetsInBarrel = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("updatedPatJetsTransientCorrectedXTag"),
    cut = cms.string("pt > 15"),
)

addModule(process.selectJetsInBarrel)

process.DeepFlavourProducer = cms.EDProducer("DeepFlavourTagInfoProducer", 
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    secondary_vertices = cms.InputTag("slimmedSecondaryVertices"),
    jets = cms.InputTag("selectJetsInBarrel"),
    shallow_tag_infos = cms.InputTag("pfDeepCSVTagInfosXTag")
)

process.nanoxProducer = cms.EDProducer("NANOXProducer",
    plugins = cms.PSet(
        globalVars = cms.PSet(
            type = cms.string("GlobalJetTagData"),
            jets = cms.InputTag("selectJetsInBarrel"),
        ),
        csv = cms.PSet(
            type = cms.string("CSVInputTagData"),
            jets = cms.InputTag("selectJetsInBarrel"),
            tagName = cms.string('pfDeepCSV')
        ),
        cpf = cms.PSet(
            type = cms.string("ChargedPFTagData"),
            jets = cms.InputTag("selectJetsInBarrel"),
            pvVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
            svVertices = cms.InputTag("slimmedSecondaryVertices"),
        ),
        npf = cms.PSet(
            type = cms.string("NeutralPFTagData"),
            jets = cms.InputTag("selectJetsInBarrel"),
            svVertices = cms.InputTag("slimmedSecondaryVertices"),
        ),
        sv = cms.PSet(
            type = cms.string("SVTagData"),
            jets = cms.InputTag("selectJetsInBarrel"),
            pvVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
            svVertices = cms.InputTag("slimmedSecondaryVertices"),
        ),
        legacyTag = cms.PSet(
            type = cms.string("LegacyTagData"),
            jets = cms.InputTag("selectJetsInBarrel"),
            pvVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
        ),
    )
)

process.nanoxFlatTable = cms.EDProducer("NANOXFlatTableProducer",
    tagData = cms.VPSet([
        cms.PSet(
            src = cms.InputTag("nanoxProducer","globalVars"),
            arrayNames = cms.vstring(["global"])
        ),
        cms.PSet(
            src = cms.InputTag("nanoxProducer","csv"),
            arrayNames = cms.vstring(["csv"])
        ),
        cms.PSet(
            src = cms.InputTag("nanoxProducer","cpf"),
            arrayNames = cms.vstring(["cpflength","cpf"])
        ),
        cms.PSet(
            src = cms.InputTag("nanoxProducer","npf"),
            arrayNames = cms.vstring(["npflength","npf"])
        ),
        cms.PSet(
            src = cms.InputTag("nanoxProducer","sv"),
            arrayNames = cms.vstring(["svlength","sv"])
        ),
        cms.PSet(
            src = cms.InputTag("nanoxProducer","legacyTag"),
            arrayNames = cms.vstring(["legacyTag"])
        ),
        
    ])
)

if not options.isData:
    process.nanoxProducer.plugins.origin = cms.PSet(
        type = cms.string("JetOriginTagData"),
        jets = cms.InputTag("selectJetsInBarrel"),
        displacedGenVertices = cms.InputTag("displacedGenVertices"),
    )
    
    process.nanoxFlatTable.tagData.append(cms.PSet(
        src = cms.InputTag("nanoxProducer","origin"),
        arrayNames = cms.vstring(["jetorigin"])
    ))
    
    if options.addLLPInfo:
        process.nanoxProducer.plugins.llpinfo = cms.PSet(
            type = cms.string("LLPInfo"),
            displacedGenVertices = cms.InputTag("displacedGenVertices"),
            LLPtype = cms.string(options.LLPtype)
        )

        process.nanoxFlatTable.tagData.append(cms.PSet(
            src = cms.InputTag("nanoxProducer","llpinfo"),
            arrayNames = cms.vstring(["llpinfo"])
        ))
    
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
 
if options.isData:
    addModule(process.nanoSequence)
else:
    addModule(process.nanoSequenceMC)
    process.load('NANOX.DisplacedVertex.GenDisplacedVertices_cff')
    addModule(process.DisplacedGenVertexSequence)

process.nanoxSequence = cms.Sequence( 
    process.DeepFlavourProducer
    +process.nanoxProducer
    +process.nanoxFlatTable
)
   

addModule(process.nanoxSequence)

if options.isData:
    process.nanoAOD_step = cms.Path(process.nanoSequence)
else:
    process.nanoAOD_step = cms.Path(process.nanoSequenceMC)
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
