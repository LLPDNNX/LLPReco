#ifndef LLPReco_DataFormats_ShallowTagInfoFeatures_h
#define LLPReco_DataFormats_ShallowTagInfoFeatures_h

namespace llpdnnx {

struct ShallowTagInfoFeatures 
{
    // jet general
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
    
    ShallowTagInfoFeatures():
        trackSumJetEtRatio(0),      // ratio of track sum transverse energy over jet energy
        trackSumJetDeltaR(0),       // pseudoangular distance between jet axis and track fourvector sum
        vertexCategory(0),          // category of secondary vertex (Reco, Pseudo, No)
        trackSip2dValAboveCharm(0), // track 2D signed impact parameter of first track lifting mass above charm
        trackSip2dSigAboveCharm(0), // track 2D signed impact parameter significance of first track lifting mass above charm
        trackSip3dValAboveCharm(0), // track 3D signed impact parameter of first track lifting mass above charm
        trackSip3dSigAboveCharm(0),  // track 3D signed impact parameter significance of first track lifting mass above charm
        // track info
        jetNTracksEtaRel(0), // tracks associated to jet for which trackEtaRel is calculated
        jetNSelectedTracks(0)
    {}

};

}

#endif
