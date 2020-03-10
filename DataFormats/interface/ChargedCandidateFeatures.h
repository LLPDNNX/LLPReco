#ifndef LLPReco_DataFormats_ChargedCandidateFeatures_h
#define LLPReco_DataFormats_ChargedCandidateFeatures_h

namespace llpdnnx {

struct ChargedCandidateFeatures
{

    int cpf_jetIdx;
    float cpf_ptrel;
    float cpf_deta;
    float cpf_dphi;

    float cpf_trackEtaRel;
    float cpf_trackPtRel;
    float cpf_trackPPar;
    float cpf_trackDeltaR;
    float cpf_trackPParRatio;
    float cpf_trackPtRatio;
    float cpf_trackSip2dVal;
    float cpf_trackSip2dSig;
    float cpf_trackSip3dVal;
    float cpf_trackSip3dSig;
    float cpf_trackJetDistVal;
    float cpf_trackJetDistSig;
    float cpf_drminsv;
    int cpf_vertex_association;
    float cpf_fromPV;
    float cpf_puppi_weight;
    float cpf_track_chi2;
    float cpf_track_quality;
    int cpf_track_numberOfValidPixelHits ;
    int cpf_track_pixelLayersWithMeasurement ;
    int cpf_track_numberOfValidStripHits ;
    int cpf_track_stripLayersWithMeasurement; 
    float cpf_relmassdrop;

    int cpf_matchedMuon;
    int cpf_matchedElectron;
    int cpf_matchedSV;
    int cpf_track_ndof;

    float cpf_dZmin;
};

}

#endif
