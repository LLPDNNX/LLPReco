function run_test()
{
    cd ~
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    export SCRAM_ARCH=slc6_amd64_gcc700 || return 1
    scramv1 project CMSSW CMSSW_10_2_18 || return 1
    cd CMSSW_10_2_18/src || return 1
    eval `scram runtime -sh` || return 1
    rsync -r /scripts/* LLPReco  || return 1
    scram b || return 1
    pip install --user -r LLPReco/requirements.txt || return 1
    cmsRun LLPReco/XTagInfoProducer/test/testXTag.py inputFiles=https://github.com/LLPDNNX/test-files/raw/master/miniaod/Moriond17_aug2018_miniAODv3_HNL.root || return 1
    cmsRun LLPReco/NANOProducer/test/produceNANO.py inputFiles=https://github.com/LLPDNNX/test-files/raw/master/miniaod/HNL_miniaod16v3.root year=2016 test=True || return 1
    python LLPReco/test/check_nan.py || return 1
    cmsRun LLPReco/NANOProducer/test/produceNANO.py inputFiles=https://github.com/LLPDNNX/test-files/raw/master/miniaod/HNL_miniaod17v2.root year=2017 test=True || return 1
    python LLPReco/test/check_nan.py || return 1
    cmsRun LLPReco/NANOProducer/test/produceNANO.py inputFiles=https://github.com/LLPDNNX/test-files/raw/master/miniaod/HNL_miniaod18.root year=2018 test=True || return 1
    python LLPReco/test/check_nan.py || return 1

}

run_test
