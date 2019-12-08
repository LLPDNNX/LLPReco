#ifndef LLPReco_DataFormats_XTagFeatures_h
#define LLPReco_DataFormats_XTagFeatures_h

#include <vector>
#include "LLPReco/DataFormats/interface/JetFeatures.h"
#include "LLPReco/DataFormats/interface/SecondaryVertexFeatures.h"
#include "LLPReco/DataFormats/interface/ShallowTagInfoFeatures.h"
#include "LLPReco/DataFormats/interface/NeutralCandidateFeatures.h"
#include "LLPReco/DataFormats/interface/ChargedCandidateFeatures.h"

namespace llpdnnx {

class XTagFeatures {

  public:
    JetFeatures jet_features;
    ShallowTagInfoFeatures tag_info_features;
    std::vector<SecondaryVertexFeatures> sv_features;
    std::vector<NeutralCandidateFeatures> npf_features;
    std::vector<ChargedCandidateFeatures> cpf_features;
    std::size_t npv; // used by deep flavour
};    


}  

#endif
