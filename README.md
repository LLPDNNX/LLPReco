# LLPReco
On-the-fly LLPDNNX tagger evaluation in miniAOD

## Checkout instructions

Using recent CMSSW release which is shipped with all the required software (especially tensorflow C++ API v1.6).
```
cmsrel CMSSW_10_2_18
cd CMSSW_10_2_18/src
cmsenv
git clone https://github.com/LLPDNNX/LLPReco.git LLPReco
scram b
```

Example configuration file
```
cmsRun LLPReco/XTagInfoProducer/test/testXTag.py 
```
