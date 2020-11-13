#ifndef LLPReco_DataFormats_ChargedCandidateFeatures_h
#define LLPReco_DataFormats_ChargedCandidateFeatures_h

#include <cmath>

namespace llpdnnx {

struct ChargedCandidateFeatures
{

    int jetIdx;
    float ptrel;
    float deta;
    float dphi;
    float deltaR;
    
    float px;
    float py;
    float pz;

    float trackEtaRel;
    float trackPtRel;
    float trackPPar;
    float trackDeltaR;
    float trackPParRatio;
    float trackPtRatio;
    float trackSip2dVal;
    float trackSip2dSig;
    float trackSip3dVal;
    float trackSip3dSig;
    float trackJetDistVal;
    float trackJetDistSig;
    float drminsv;
    int vertex_association;
    float fromPV;
    float puppi_weight;
    float track_chi2;
    float track_quality;
    int track_numberOfValidPixelHits;
    int track_pixelLayersWithMeasurement;
    int track_numberOfValidStripHits;
    int track_stripLayersWithMeasurement; 
    float relmassdrop;

    int matchedMuon;
    int matchedElectron;
    int matchedSV;
    int track_ndof;

    float dZmin;
    
    ChargedCandidateFeatures():
        jetIdx(-1),
        ptrel(-1),
        deta(1),
        dphi(1),
        deltaR(1),

        px(0),
        py(0),
        pz(0),

        trackEtaRel(0),
        trackPtRel(0),
        trackPPar(0),
        trackDeltaR(0),
        trackPParRatio(0),
        trackPtRatio(0),
        trackSip2dVal(0),
        trackSip2dSig(0),
        trackSip3dVal(0),
        trackSip3dSig(0),
        trackJetDistVal(0),
        trackJetDistSig(0),
        drminsv(0),
        vertex_association(0),
        fromPV(0),
        puppi_weight(0),
        track_chi2(0),
        track_quality(0),
        track_numberOfValidPixelHits(0),
        track_pixelLayersWithMeasurement(0),
        track_numberOfValidStripHits(0),
        track_stripLayersWithMeasurement(0),
        relmassdrop(0),

        matchedMuon(0),
        matchedElectron(0),
        matchedSV(0),
        track_ndof(0),

        dZmin(100)
    {}
    
    bool operator<(const ChargedCandidateFeatures& other) const
    {
        if (trackSip2dSig>0 and other.trackSip2dSig>0)
        {
            return std::fabs(trackSip2dSig)>std::fabs(other.trackSip2dSig); //sort decreasing
        }
        else if (trackSip2dSig<0 and other.trackSip2dSig>0)
        {
            return false;
        }
        else if (trackSip2dSig>0 and other.trackSip2dSig<0)
        {
            return true;
        }
        else if (std::fabs(drminsv-other.drminsv)>std::numeric_limits<float>::epsilon())
        {
            return drminsv<other.drminsv; //sort increasing
        }
        else
        {
            return ptrel>other.ptrel;  //sort decreasing
        }
        
        return false;
    }
};

}

#endif
