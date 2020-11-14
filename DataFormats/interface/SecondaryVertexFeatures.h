#ifndef LLPReco_DataFormats_SecondaryVertexFeatures_h
#define LLPReco_DataFormats_SecondaryVertexFeatures_h

#include <cmath>

namespace llpdnnx {

struct SecondaryVertexFeatures
{

    int jetIdx;
    float ptrel;
    float deta;
    float dphi;
    float deltaR;
    float mass;
    int ntracks;
    float chi2;
    int ndof;
    float dxy;
    float dxysig;
    float d3d;
    float d3dsig;
    float costhetasvpv;
    float enratio;
    
    bool operator<(const SecondaryVertexFeatures& other) const
    {
        if (dxysig>0 and other.dxysig>0)
        {
            if (std::fabs(dxysig-other.dxysig)>std::numeric_limits<float>::epsilon())
            {
                return dxysig>other.dxysig;
            }
        }
        return ptrel>other.ptrel; //sort decreasing
    }

};

}

#endif
