#ifndef LLPReco_DataFormats_ShallowTagInfoFeatures_h
#define LLPReco_DataFormats_ShallowTagInfoFeatures_h

namespace llpdnnx {

struct ShallowTagInfoFeatures 
{
    // jet general
    int jetIdx;
    float trackSumJetEtRatio;      // ratio of track sum transverse energy over jet energy
    float trackSumJetDeltaR;       // pseudoangular distance between jet axis and track fourvector sum
    int vertexCategory;          // category of secondary vertex (Reco, Pseudo, No)
    float trackSip2dValAboveCharm; // track 2D signed impact parameter of first track lifting mass above charm
    float trackSip2dSigAboveCharm; // track 2D signed impact parameter significance of first track lifting mass above charm
    float trackSip3dValAboveCharm; // track 3D signed impact parameter of first track lifting mass above charm
    float trackSip3dSigAboveCharm; // track 3D signed impact parameter significance of first track lifting mass above charm
    // track info
    int jetNTracksEtaRel; // tracks associated to jet for which trackEtaRel is calculated
    int jetNSelectedTracks;

};

}

#endif
