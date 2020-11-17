#ifndef LLPReco_DataFormats_NeutralCandidateFeatures_h
#define LLPReco_DataFormats_NeutralCandidateFeatures_h

#include <cmath>

namespace llpdnnx {

struct NeutralCandidateFeatures 
{
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
    
    NeutralCandidateFeatures():
        ptrel(0),
        deta(0),
        dphi(0),
        deltaR(0),

        px(0),
        py(0),
        pz(0),

        isGamma(0),
        hcal_fraction(0),
        drminsv(0),
        puppi_weight(0),
        relmassdrop(0)
    {}
    
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
