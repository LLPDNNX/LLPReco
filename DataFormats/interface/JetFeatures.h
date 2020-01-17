#ifndef LLPReco_DataFormats_JetFeatures_h
#define LLPReco_DataFormats_JetFeatures_h

namespace llpdnnx {

class JetFeatures {

  public:

    int jetIdx;
    float pt;
    float eta;
    float mass;
    float energy;

    float tau1;
    float tau2;
    float tau3;
    
    float massDropMassAK;
    float massDropMassCA;
    float softDropMassAK;
    float softDropMassCA;
    
    float thrust;
    float sphericity;
    float circularity;
    float isotropy;
    float eventShapeC;
    float eventShapeD;
};

}

#endif 
