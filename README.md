# LLPReco
* On-the-fly LLPDNNX tagger evaluation in miniAOD
* Production of extended nanoAOD ntuples for analysis and training

![build tests](https://travis-ci.org/LLPDNNX/LLPReco.svg?branch=master)


## Checkout instructions

Using recent CMSSW release which is shipped with all the required software (especially tensorflow C++ API v1.6).
```
cmsrel CMSSW_10_2_18
cd CMSSW_10_2_18/src
cmsenv
git clone https://github.com/LLPDNNX/LLPReco.git LLPReco
scram b
```

Configuration file for full Run 2 extended nanoAOD ntuple production
```
cmsRun LLPReco/NANOProducer/test/produceNANO.py 
```

Example configuration file to run the tagger evaluation
```
cmsRun LLPReco/XTagInfoProducer/test/testXTag.py 
```

The output can be analysed in a similar way to b-tagging, as described here <https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookBTagging>

For example, ```recoJetedmRefToBaseProdTofloatsAssociationVector_pfXTags_1000_USER.obj.data_``` will correspond to the tagger evaluated at 1000 mm. Note that tagger output might be meaningless for jets with |η| > 2.4.
