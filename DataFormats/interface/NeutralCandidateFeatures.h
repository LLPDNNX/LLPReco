#ifndef LLPReco_DataFormats_NeutralCandidateFeatures_h
#define LLPReco_DataFormats_NeutralCandidateFeatures_h

#include <cmath>

namespace llpdnnx {

struct NeutralCandidateFeatures 
{
    int jetIdx;
    float ptrel;
    float deta;
    float dphi;
    float deltaR;
    
    float px;
    float py;
    float pz;
    
    int isGamma;
    float hcal_fraction;
    float drminsv;
    float puppi_weight;
    float relmassdrop;
    
    bool operator<(const NeutralCandidateFeatures& other) const
    {
        if (std::fabs(drminsv-other.drminsv)>std::numeric_limits<float>::epsilon())
        {
            return drminsv<other.drminsv; //sort increasing
        }
        else
        {
            return ptrel>other.ptrel; //sort decreasing
        }
        return false;
    }
};

}

#endif
