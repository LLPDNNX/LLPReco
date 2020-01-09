#ifndef LLPReco_DataFormats_XTagInfo_h
#define LLPReco_DataFormats_XTagInfo_h

#include "LLPReco/DataFormats/interface/XTagFeatures.h"
#include "DataFormats/BTauReco/interface/FeaturesTagInfo.h"

namespace reco {

    typedef  FeaturesTagInfo<llpdnnx::XTagFeatures> XTagInfo;

    DECLARE_EDM_REFS( XTagInfo )

}

#endif

