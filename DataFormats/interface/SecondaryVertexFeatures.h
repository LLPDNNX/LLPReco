#ifndef LLPReco_DataFormats_SecondaryVertexFeatures_h
#define LLPReco_DataFormats_SecondaryVertexFeatures_h

namespace llpdnnx {

class SecondaryVertexFeatures {

  public:

    float sv_pt;
    float sv_deltaR;
    float sv_mass;
    float sv_ntracks;
    float sv_chi2;
    float sv_ndof;
    float sv_dxy;
    float sv_dxysig;
    float sv_d3d;
    float sv_d3dsig;
    float sv_costhetasvpv;
    float sv_enratio;

};

}

#endif 
