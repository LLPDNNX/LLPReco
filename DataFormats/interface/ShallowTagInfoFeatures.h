#ifndef LLPReco_DataFormats_ShallowTagInfoFeatures_h
#define LLPReco_DataFormats_ShallowTagInfoFeatures_h

namespace llpdnnx {

struct ShallowTagInfoFeatures 
{
    // jet general
    int csv_jetIdx;
    float csv_trackSumJetEtRatio;      // ratio of track sum transverse energy over jet energy
    float csv_trackSumJetDeltaR;       // pseudoangular distance between jet axis and track fourvector sum
    float csv_vertexCategory;          // category of secondary vertex (Reco, Pseudo, No)
    float csv_trackSip2dValAboveCharm; // track 2D signed impact parameter of first track lifting mass above charm
    float csv_trackSip2dSigAboveCharm; // track 2D signed impact parameter significance of first track lifting mass above charm
    float csv_trackSip3dValAboveCharm; // track 3D signed impact parameter of first track lifting mass above charm
    float csv_trackSip3dSigAboveCharm; // track 3D signed impact parameter significance of first track lifting mass above charm
    // track info
    int csv_jetNTracksEtaRel; // tracks associated to jet for which trackEtaRel is calculated
    int csv_jetNSelectedTracks;

};

}

#endif
