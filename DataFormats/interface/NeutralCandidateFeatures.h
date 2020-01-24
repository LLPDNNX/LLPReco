#ifndef LLPReco_DataFormats_NeutralCandidateFeatures_h
#define LLPReco_DataFormats_NeutralCandidateFeatures_h

namespace llpdnnx {

struct NeutralCandidateFeatures 
{
    int npf_jetIdx;
    float npf_ptrel;
    float npf_deta;
    float npf_dphi;
    float npf_deltaR;
    float npf_isGamma;
    float npf_hcal_fraction;
    float npf_drminsv;
    float npf_puppi_weight;
    float npf_relmassdrop;
};

}

#endif
