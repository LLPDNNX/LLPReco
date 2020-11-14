#ifndef LLPReco_DataFormats_MuonCandidateFeatures_h
#define LLPReco_DataFormats_MuonCandidateFeatures_h

namespace llpdnnx {

struct MuonCandidateFeatures {

    int jetIdx;
    int isGlobal;
    int isTight;
    int isMedium;
    int isLoose;
    int isStandAlone;

    float ptrel;
    float deta;
    float dphi;
    float px;
    float py;
    float pz;
    float charge;
    float energy;
    float et;
    float deltaR;
    int numberOfMatchedStations;

    float IP2d;
    float IP2dSig;
    float IP3d;
    float IP3dSig;

    float EtaRel;
    float dxy;
    float dxyError;
    float dxySig;
    float dz;
    float dzError;
    float dzSig;
    int numberOfValidPixelHits;
    int numberOfpixelLayersWithMeasurement;
    int numberOfstripLayersWithMeasurement; //that does not help. needs to be discussed.

    float chi2;
    int ndof;

    float caloIso;
    float ecalIso;
    float hcalIso;

    float sumPfChHadronPt;
    float sumPfNeuHadronEt;
    float Pfpileup;
    float sumPfPhotonEt;

    float sumPfChHadronPt03;
    float sumPfNeuHadronEt03;
    float Pfpileup03;
    float sumPfPhotonEt03;


    float timeAtIpInOut;
    float timeAtIpInOutErr;
    float timeAtIpOutIn;
    
    
    MuonCandidateFeatures():
        jetIdx(-1),
        isGlobal(-1),
        isTight(-1),
        isMedium(-1),
        isLoose(-1),
        isStandAlone(-1),

        ptrel(0),
        deta(0),
        dphi(0),
        px(0),
        py(0),
        pz(0),
        charge(0),
        energy(0),
        et(0),
        deltaR(0),
        numberOfMatchedStations(0),

        IP2d(0),
        IP2dSig(0),
        IP3d(0),
        IP3dSig(0),

        EtaRel(0),
        dxy(0),
        dxyError(0),
        dxySig(0),
        dz(0),
        dzError(0),
        dzSig(0),
        numberOfValidPixelHits(0),
        numberOfpixelLayersWithMeasurement(0),
        numberOfstripLayersWithMeasurement(0), //that does not help. needs to be discussed.

        chi2(0),
        ndof(0),

        caloIso(0),
        ecalIso(0),
        hcalIso(0),

        sumPfChHadronPt(0),
        sumPfNeuHadronEt(0),
        Pfpileup(0),
        sumPfPhotonEt(0),

        sumPfChHadronPt03(0),
        sumPfNeuHadronEt03(0),
        Pfpileup03(0),
        sumPfPhotonEt03(0),


        timeAtIpInOut(0),
        timeAtIpInOutErr(0),
        timeAtIpOutIn(0)
    {}

    bool operator<(const MuonCandidateFeatures& other) const
    {
        if (IP2dSig>0 and other.IP2dSig>0)
        {
            if (std::fabs(IP2dSig-other.IP2dSig)>std::numeric_limits<float>::epsilon())
            {
                return std::fabs(IP2dSig)>std::fabs(other.IP2dSig); //sort decreasing
            }
        }
        return ptrel>other.ptrel; //sort decreasing
    }

};


}

#endif
