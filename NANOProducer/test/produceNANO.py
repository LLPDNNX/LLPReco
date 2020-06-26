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
    'addSignalLHE',
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "adds LHE weights of signal samples"
)

options.register(
    'addLLPInfo',
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "add LLP Info"
)

options.register(
    'year',
    '2016',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "add year file"
)
options.parseArguments() 

if options.year == '2016':
    process = cms.Process('NANO',eras.Run2_2016,eras.run2_nanoAOD_94X2016)
elif options.year == '2017':
    process = cms.Process('NANO',eras.Run2_2017,eras.run2_nanoAOD_94XMiniAODv2)
elif options.year == '2018' or options.year == '2018D':
    process = cms.Process('NANO',eras.Run2_2018,eras.run2_nanoAOD_102Xv1)
else:
    process = cms.Process('NANO',eras.Run2_2016,eras.run2_nanoAOD_94X2016)

print "Selected year: ",options.year

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


if options.isData:
    process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
    dataTier = cms.untracked.string('NANOAOD')
else:
    process.load('SimGeneral.MixingModule.mixNoPU_cfi')
    process.load('Configuration.StandardSequences.MagneticField_cff')
    dataTier = cms.untracked.string('NANOAODSIM')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.options = cms.untracked.PSet()

files = {
    'test': {
        "mc": "/store/user/kjpena/miniAODv3_08Feb2020/GluGluH_HToSSTobbbb_MH-125_MS-25_ctauS-500_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_MINIAODSIM/200209_212810/0000/output_1.root",
        #"mc": "https://github.com/LLPDNNX/test-files/raw/master/miniaod/Moriond17_aug2018_miniAODv3_HNL.root",
        },
    '2016': {
        "mc": "root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mkomm/HNL/miniaod16v3_200625/HNL_dirac_all_ctau1p0e-02_massHNL8p0_Vall2p996e-02/miniaod16v3_200625/200625_171726/0000/HNL2016_9.root",
        #"mc": "root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mkomm/HNL/miniaod16v3_200517/HNL_dirac_all_ctau1p0e00_massHNL6p0_Vall6p496e-03/miniaod16v3_200517/200517_004822/0000/HNL2016_140.root", #"root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-8_V-0.004472135955_tau_Dirac_massiveAndCKM_LO/heavyNeutrino_1.root",
        #"mc": "root://maite.iihe.ac.be///store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-10_V-0.00112249721603_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_76.root",
        "data": "/store/data/Run2016H/SingleMuon/MINIAOD/17Jul2018-v1/00000/16924A85-4D8C-E811-A51C-A4BF01013F29.root",
    },
    '2017': {
        "mc": "root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-8_V-0.00214242852856_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_10.root",
        "data": "/store/data/Run2017E/SingleMuon/MINIAOD/31Mar2018-v1/00000/A6325FCE-1C39-E811-BB22-0CC47A745298.root"
    },
    '2018': {
        "mc": "root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Autumn18/displaced/HeavyNeutrino_lljj_M-8_V-0.00214242852856_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_10.root",
        "data": "/store/data/Run2018B/SingleMuon/MINIAOD/17Sep2018-v1/60000/FF47BB90-FC1A-CC44-A635-2B8B8C64AA39.root"
    },
    '2018D': {
        "data": "/store/data/Run2018B/SingleMuon/MINIAOD/17Sep2018-v1/60000/FF47BB90-FC1A-CC44-A635-2B8B8C64AA39.root"
    }
}

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(files[options.year]['data'] if options.isData else files[options.year]['mc'])

)
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('test102X nevts:10000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)


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
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('llpnanoAOD_step') #only events passing this path will be saved
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
        
        'drop *_saJetTable_*_*',
        'drop *_FatJetTable_*_*',
        'drop *_saTable_*_*',
        
        'drop *_rivetMetTable_*_*',
        'drop *_rivetProducerHTXS_*_*',
        
        #'drop *_rivetMetTable_*_*',
    )
)

## Output file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.OUT = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    outputCommands = cms.untracked.vstring(['keep *'])
)

if options.year == "test":
    options.year = "2016"

if options.isData:
    if options.year == '2016':
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v12', '')
    if options.year == '2017':
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v12', '')
    if options.year == '2018':
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v12', '')
    if options.year == '2018D':
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Prompt_v15', '')
    jetCorrectionsAK4PFchs = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'], 'None')
else:
    if options.year == '2016':
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_mcRun2_asymptotic_v7', '')
    if options.year == '2017':
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_mc2017_realistic_v7', '')
    if options.year == '2018' or options.year == '2018D':
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v20', '')
    jetCorrectionsAK4PFchs = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')


from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
    process,
    labelName = "XTag",
    jetSource = cms.InputTag('updatedJets'),
    jetCorrections = jetCorrectionsAK4PFchs,
    pfCandidates = cms.InputTag('packedPFCandidates'),
    pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
    #svSource = cms.InputTag('adaptedSlimmedSecondaryVertices'), 
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    muSource = cms.InputTag('slimmedMuons'),
    elSource = cms.InputTag('slimmedElectrons'),
    btagInfos = [
        'pfImpactParameterTagInfos',
        'pfInclusiveSecondaryVertexFinderTagInfos',
        'pfDeepCSVTagInfos',
    ],
    btagDiscriminators = ['pfCombinedInclusiveSecondaryVertexV2BJetTags'],
    explicitJTA = False,
)

process.pfXTagInfos = cms.EDProducer("XTagInfoProducer",
    jets = cms.InputTag("updatedJets"),
    muonSrc  = cms.InputTag("slimmedMuons"),
    electronSrc = cms.InputTag("slimmedElectrons"),
    shallow_tag_infos = cms.InputTag('pfDeepCSVTagInfosXTag'),
    vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
    #secondary_vertices = cms.InputTag("adaptedSlimmedSecondaryVertices")
    secondary_vertices = cms.InputTag("slimmedSecondaryVertices")
)

process.nanoTable = cms.EDProducer("NANOProducer",
    srcJets = cms.InputTag("finalJets"),
    srcTags = cms.InputTag("pfXTagInfos"),
)

process.nanoGenTable = cms.EDProducer("NANOGenProducer",
    srcJets = cms.InputTag("finalJets"),
    srcLabels = cms.InputTag("llpLabels"),
    srcTags = cms.InputTag("pfXTagInfos")
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
    srcJets = cms.InputTag("updatedJets"),
    srcDecayInfo = cms.InputTag("llpGenDecayInfo"),
    jetAlgorithm = cms.string("AntiKt"),
    rParam = cms.double(0.4),
    ghostRescaling = cms.double(1e-18),
    relPtTolerance = cms.double(1e-3)
)

process.llpLabels = cms.EDProducer(
    "LLPLabelProducer",
    srcVertices = cms.InputTag("displacedGenVertices"),
    srcJets = cms.InputTag("updatedJets"),
    srcFlavourInfo = cms.InputTag("llpFlavour"),
    tauPtThreshold = cms.double(1.),
    quarkPtThreshold = cms.double(1.),
    bPtThreshold = cms.double(1.),
    muonPtThreshold = cms.double(1.),
    electronPtThreshold = cms.double(1.),
)



process.lheWeightsTable = cms.EDProducer(
    "LHEWeightsProducer",
    lheInfo = cms.VInputTag(cms.InputTag("externalLHEProducer"), cms.InputTag("source")),
    weightGroups = cms.PSet()
)

#coupling reweighting
process.lheWeightsTable.weightGroups.coupling = cms.vstring()
for i in range(1,68):
    process.lheWeightsTable.weightGroups.coupling.append("rwgt_%i"%(i))
    
#PDF NNPDF3.1 NNLO hessian
process.lheWeightsTable.weightGroups.nnpdfhessian = cms.vstring()
for i in range(1048,1151):
    process.lheWeightsTable.weightGroups.nnpdfhessian.append("%i"%(i))
    
#PDF NNPDF3.1 NNLO replicas
process.lheWeightsTable.weightGroups.nnpdfreplica = cms.vstring()
for i in range(1151,1252):
    process.lheWeightsTable.weightGroups.nnpdfreplica.append("%i"%(i))
    
#scale weights
for scaleSet in [
    ['murNominal_mufNominal',range(1001,1006)],
    ['murUp_mufNominal',range(1006,1011)],
    ['murDown_mufNominal',range(1011,1016)],
    ['murNominal_mufUp',range(1016,1021)],
    ['murUp_mufUp',range(1021,1026)],
    ['murDown_mufUp',range(1026,1031)],
    ['murNominal_mufDown',range(1031,1036)],
    ['murUp_mufDown',range(1036,1041)],
    ['murDown_mufDown',range(1041,1046)],
    ['emission',range(1046,1048)],
]:
    
    setattr(process.lheWeightsTable.weightGroups,scaleSet[0],cms.vstring())
    for i in scaleSet[1]:
        getattr(process.lheWeightsTable.weightGroups,scaleSet[0]).append("%i"%(i))
        
        

#process.load('RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff')
#process.load('LLPReco.NANOProducer.adaptedSV_cff')


process.selectedMuonsForFilter = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string("pt>25 && isGlobalMuon()")
)
process.selectedMuonsMinFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("selectedMuonsForFilter"),
    minNumber = cms.uint32(1)
)
    
process.muonFilterSequence = cms.Sequence(
    process.selectedMuonsForFilter+process.selectedMuonsMinFilter
)

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeData,nanoAOD_customizeMC

#call to customisation function nanoAOD_customizeMC imported from PhysicsTools.NanoAOD.nano_cff
if options.isData:
    process = nanoAOD_customizeData(process)
else:
    process = nanoAOD_customizeMC(process)
    


if options.isData:
    process.llpnanoAOD_step = cms.Path(
        process.muonFilterSequence+
        process.nanoSequence+
        #process.adaptedVertexing+
        process.pfXTagInfos+
        process.nanoTable
    )
else:
    process.llpnanoAOD_step = cms.Path(
        process.nanoSequenceMC+
        #process.adaptedVertexing+
        process.pfXTagInfos+
        process.displacedGenVertexSequence+
        process.llpGenDecayInfo+
        process.llpFlavour+
        process.llpLabels+
        process.nanoTable+
        process.nanoGenTable
    )
    
    if options.addSignalLHE:
        process.llpnanoAOD_step += process.lheWeightsTable
    
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)



process.schedule = cms.Schedule(process.llpnanoAOD_step,process.endjob_step,process.NANOAODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


modulesToRemove = [
    'jetCorrFactorsAK8',
    'updatedJetsAK8',
    'tightJetIdAK8',
    'looseJetIdAK8',
    'tightJetIdLepVetoAK8',
    'updatedJetsAK8WithUserData',
    'chsForSATkJets',
    'softActivityJets',
    'softActivityJets2',
    'softActivityJets5',
    'softActivityJets10',
    'finalJetsAK8',
    'fatJetTable',
    'subJetTable',
    'saJetTable',
    'saTable',


    "genJetAK8Table",
    "genJetAK8FlavourAssociation",
    "genJetAK8FlavourTable",
   
    "HTXSCategoryTable",
    "rivetProducerHTXS",
    "genSubJetAK8Table",
    
    "l1bits",
]

#override final jets

#process.finalJets.addBTagInfo=cms.bool(True)
#process.finalJets.addDiscriminators = cms.bool(True)
#process.finalJets.addTagInfos=cms.bool(True)

#remove unneeded modules
for moduleName in modulesToRemove:
    if hasattr(process,moduleName):
        print "removing module: ",moduleName
        if options.isData:
            process.nanoSequence.remove(getattr(process,moduleName))
        else:
            process.nanoSequenceMC.remove(getattr(process,moduleName))
    else:
        print "module for removal not found: ",moduleName

#override final photons (required by object linker) so that ID evaluation is not needed
#process.finalPhotons.cut = cms.string("pt > 5")
#process.finalPhotons.src = cms.InputTag("slimmedPhotons")

'''
process.MINIAODoutput = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('output.root'),
    outputCommands = process.NANOAODSIMoutput.outputCommands,
    dropMetaData = cms.untracked.string('ALL'),
)
'''
    
process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
#print process.dumpPython()
# End adding early deletion

#process.endpath= cms.EndPath(process.OUT)
