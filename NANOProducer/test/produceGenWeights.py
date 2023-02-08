import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
#from PhysicsTools.NanoAOD.common_cff import *
from Configuration.StandardSequences.Eras import eras
from RecoBTag.Configuration.RecoBTag_cff import *
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process('GEN',eras.Run2_2016)

import random

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic50ns13TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')


options = VarParsing ('analysis')

options.register(
    'gridpack',
    '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "is data"
)

options.register(
    'output',
    '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "is data"
)
options.parseArguments() 


process.RandomNumberGeneratorService.generator.initialSeed = random.randint(1,10000000)
process.RandomNumberGeneratorService.externalLHEProducer.initialSeed = random.randint(1,10000000)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(8000)
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source("EmptySource")


#Output definition
process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    saveProvenance = cms.untracked.bool(True),
    fakeNameForCrab = cms.untracked.bool(True),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
            ['nanoAOD_step']
        ) #only events passing this path will be saved
    ),
    fileName = cms.untracked.string(options.output+'.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep nanoaodFlatTable_*Table_*_*',
        'keep edmTriggerResults_*_*_*',
        'keep nanoaodMergeableCounterTable_*Table_*_*',
        'keep nanoaodUniqueString_nanoMetadata_*_*',
        
        'drop *_caloMetTable_*_*',
        'drop *_saJetTable_*_*',
        'drop *_saTable_*_*',
        'drop *_fatJetTable_*_*',
        'drop *_fatJetMCTable_*_*',
        'drop *_subJetTable_*_*',
        'drop *_subjetMCTable_*_*',
        'drop *_genJetAK8FlavourTable_*_*',
        'drop *_genJetAK8Table_*_*',
        'drop *_genSubJetAK8Table_*_*',
        'drop *_genVisTauTable_*_*',
        'drop *_tkMetTable_*_*',
        'drop *_puppiMetTable_*_*',
        'drop *_ttbarCategoryTable_*_*',
        
        
        'drop *_rivetMetTable_*_*',
        'drop *_rivetProducerHTXS_*_*',
        
        #'drop *_rivetMetTable_*_*',
    )
)


process.NANOAODSIMoutputTau0Filtered = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    saveProvenance = cms.untracked.bool(True),
    fakeNameForCrab = cms.untracked.bool(True),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
            ['tau0Filter_step']
        ) #only events passing this path will be saved
    ),
    fileName = cms.untracked.string(options.output+'_tau0.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep nanoaodFlatTable_*Table_*_*',
        'keep edmTriggerResults_*_*_*',
        'keep nanoaodMergeableCounterTable_*Table_*_*',
        'keep nanoaodUniqueString_nanoMetadata_*_*',
        
        'drop *_caloMetTable_*_*',
        'drop *_saJetTable_*_*',
        'drop *_saTable_*_*',
        'drop *_fatJetTable_*_*',
        'drop *_fatJetMCTable_*_*',
        'drop *_subJetTable_*_*',
        'drop *_subjetMCTable_*_*',
        'drop *_genJetAK8FlavourTable_*_*',
        'drop *_genJetAK8Table_*_*',
        'drop *_genSubJetAK8Table_*_*',
        'drop *_genVisTauTable_*_*',
        'drop *_tkMetTable_*_*',
        'drop *_puppiMetTable_*_*',
        'drop *_ttbarCategoryTable_*_*',
        
        
        'drop *_rivetMetTable_*_*',
        'drop *_rivetProducerHTXS_*_*',
        
        #'drop *_rivetMetTable_*_*',
    )
)

process.NANOAODSIMoutputTau1Filtered = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    saveProvenance = cms.untracked.bool(True),
    fakeNameForCrab = cms.untracked.bool(True),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
            ['tau1Filter_step']
        ) #only events passing this path will be saved
    ),
    fileName = cms.untracked.string(options.output+'_tau1.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep nanoaodFlatTable_*Table_*_*',
        'keep edmTriggerResults_*_*_*',
        'keep nanoaodMergeableCounterTable_*Table_*_*',
        'keep nanoaodUniqueString_nanoMetadata_*_*',
        
        'drop *_caloMetTable_*_*',
        'drop *_saJetTable_*_*',
        'drop *_saTable_*_*',
        'drop *_fatJetTable_*_*',
        'drop *_fatJetMCTable_*_*',
        'drop *_subJetTable_*_*',
        'drop *_subjetMCTable_*_*',
        'drop *_genJetAK8FlavourTable_*_*',
        'drop *_genJetAK8Table_*_*',
        'drop *_genSubJetAK8Table_*_*',
        'drop *_genVisTauTable_*_*',
        'drop *_tkMetTable_*_*',
        'drop *_puppiMetTable_*_*',
        'drop *_ttbarCategoryTable_*_*',
        
        
        'drop *_rivetMetTable_*_*',
        'drop *_rivetProducerHTXS_*_*',
        
        #'drop *_rivetMetTable_*_*',
    )
)

process.NANOAODSIMoutputTau2Filtered = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    saveProvenance = cms.untracked.bool(True),
    fakeNameForCrab = cms.untracked.bool(True),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
            ['tau2Filter_step']
        ) #only events passing this path will be saved
    ),
    fileName = cms.untracked.string(options.output+'_tau2.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep nanoaodFlatTable_*Table_*_*',
        'keep edmTriggerResults_*_*_*',
        'keep nanoaodMergeableCounterTable_*Table_*_*',
        'keep nanoaodUniqueString_nanoMetadata_*_*',
        
        'drop *_caloMetTable_*_*',
        'drop *_saJetTable_*_*',
        'drop *_saTable_*_*',
        'drop *_fatJetTable_*_*',
        'drop *_fatJetMCTable_*_*',
        'drop *_subJetTable_*_*',
        'drop *_subjetMCTable_*_*',
        'drop *_genJetAK8FlavourTable_*_*',
        'drop *_genJetAK8Table_*_*',
        'drop *_genSubJetAK8Table_*_*',
        'drop *_genVisTauTable_*_*',
        'drop *_tkMetTable_*_*',
        'drop *_puppiMetTable_*_*',
        'drop *_ttbarCategoryTable_*_*',
        
        
        'drop *_rivetMetTable_*_*',
        'drop *_rivetProducerHTXS_*_*',
        
        #'drop *_rivetMetTable_*_*',
    )
)

from PhysicsTools.NanoAOD.common_cff import *



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
        



   
# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeData,nanoAOD_customizeMC,ExtVar
process = nanoAOD_customizeMC(process)
    
process.nanoAOD_step = cms.Path(process.lheWeightsTable)
    
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput*process.NANOAODSIMoutputTau0Filtered*process.NANOAODSIMoutputTau1Filtered*process.NANOAODSIMoutputTau2Filtered)


process.generator = cms.EDFilter("Pythia8HadronizerFilter",
    PythiaParameters = cms.PSet(
        parameterSets = cms.vstring('pythia8CommonSettings', 
            'pythia8CUEP8M1Settings'),
        pythia8CUEP8M1Settings = cms.vstring('Tune:pp 14', 
            'Tune:ee 7', 
            'MultipartonInteractions:pT0Ref=2.4024', 
            'MultipartonInteractions:ecmPow=0.25208', 
            'MultipartonInteractions:expPow=1.6'),
        pythia8CommonSettings = cms.vstring('Tune:preferLHAPDF = 2', 
            'Main:timesAllowErrors = 10000', 
            'Check:epTolErr = 0.01', 
            'Beams:setProductionScalesFromLHEF = off', 
            'SLHA:keepSM = on', 
            'SLHA:minMassSM = 1000.', 
            'ParticleDecays:limitTau0 = on', 
            'ParticleDecays:tau0Max = 10', 
            'ParticleDecays:allowPhotonRadiation = on')
    ),
    comEnergy = cms.double(13000.0),
    filterEfficiency = cms.untracked.double(1.0),
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    pythiaPylistVerbosity = cms.untracked.int32(1)
)


process.externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    args = cms.vstring(options.gridpack),
    nEvents = cms.untracked.uint32(8000),
    numberOfParameters = cms.uint32(1),
    outputFile = cms.string('cmsgrid_final.lhe'),
    scriptName = cms.FileInPath('LLPReco/NANOProducer/data/run_generic_tarball_xrootd.sh')
)
'''
process.leptonFilter = cms.EDFilter("MCSingleParticleFilter",
    MinPt = cms.untracked.vdouble(20.0, 20.0, 20.0, 20.0),
    ParticleID = cms.untracked.vint32(11, -11, 13, -13)
)
'''

process.hnlDecayInfoTable = cms.EDProducer("HNLDecayInfoProducer",
    srcGenParticles = cms.InputTag("genParticles")
)

process.tau0Filter = cms.EDFilter("GenLeptonPairFilter",
    srcGenParticles = cms.InputTag("genParticles"),
    srcGenJets = cms.InputTag("ak4GenJetsNoNu"),
    nTaus=cms.int32(0)
)

process.tau1Filter = cms.EDFilter("GenLeptonPairFilter",
    srcGenParticles = cms.InputTag("genParticles"),
    srcGenJets = cms.InputTag("ak4GenJetsNoNu"),
    nTaus=cms.int32(1)
)

process.tau2Filter = cms.EDFilter("GenLeptonPairFilter",
    srcGenParticles = cms.InputTag("genParticles"),
    srcGenJets = cms.InputTag("ak4GenJetsNoNu"),
    nTaus=cms.int32(2)
)

process.lhe_step = cms.Path(process.externalLHEProducer)
process.generation_step = cms.Path(process.generator*process.pgen*process.hnlDecayInfoTable)
process.tau0Filter_step = cms.Path(process.tau0Filter)
process.tau1Filter_step = cms.Path(process.tau1Filter)
process.tau2Filter_step = cms.Path(process.tau2Filter)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(
    process.lhe_step,process.generation_step,
    process.tau0Filter_step,process.tau1Filter_step,process.tau2Filter_step,
    process.nanoAOD_step,process.genfiltersummary_step,process.endjob_step,
    process.NANOAODSIMoutput_step
)

process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)


