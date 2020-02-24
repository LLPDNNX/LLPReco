#ifndef LLPReco_DataFormats_LLPLabelInfo_h
#define LLPReco_DataFormats_LLPLabelInfo_h

#include "LLPReco/DataFormats/interface/LLPLabel.h"
#include "DataFormats/BTauReco/interface/FeaturesTagInfo.h"

namespace reco {

typedef FeaturesTagInfo<llpdnnx::LLPLabel> LLPLabelInfo;

DECLARE_EDM_REFS(LLPLabelInfo)

}

#endif
