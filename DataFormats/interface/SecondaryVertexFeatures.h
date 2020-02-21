#ifndef LLPReco_DataFormats_SecondaryVertexFeatures_h
#define LLPReco_DataFormats_SecondaryVertexFeatures_h

namespace llpdnnx {

struct SecondaryVertexFeatures
{

    int sv_jetIdx;
    float sv_ptrel;
    float sv_deta;
    float sv_dphi;
    float sv_deltaR;
    float sv_mass;
    int sv_ntracks;
    float sv_chi2;
    int sv_ndof;
    float sv_dxy;
    float sv_dxysig;
    float sv_d3d;
    float sv_d3dsig;
    float sv_costhetasvpv;
    float sv_enratio;

};

}

#endif
