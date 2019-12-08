#ifndef LLPReco_DataFormats_ChargedCandidateFeatures_h
#define LLPReco_DataFormats_ChargedCandidateFeatures_h

namespace llpdnnx {

class ChargedCandidateFeatures {

  public:

    float cpf_trackEtaRel;
    float cpf_trackPtRel;
    float cpf_trackPPar;
    float cpf_trackDeltaR;
    float cpf_trackPParRatio;
    //float cpf_trackPtRatio;
    float cpf_trackSip2dVal;
    float cpf_trackSip2dSig;
    float cpf_trackSip3dVal;
    float cpf_trackSip3dSig;
    float cpf_trackJetDistVal;
    float cpf_trackJetDistSig;
    float cpf_ptrel;
    float cpf_drminsv;
    float cpf_vertex_association;
    float cpf_fromPV;
    float cpf_puppi_weight;
    float cpf_track_chi2;
    float cpf_track_quality;
    float cpf_track_ndof;
};

}

#endif
