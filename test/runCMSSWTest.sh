function run_test()
{
    cd ~ 
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    export SCRAM_ARCH=slc6_amd64_gcc700 || return 1
    scramv1 project CMSSW CMSSW_10_2_18 || return 1
    cd CMSSW_10_2_18/src || return 1
    eval `scramv1 runtime -sh` || return 1
    scp -r /scripts/ . || return 1
    scram b || return 1
}

run_test
