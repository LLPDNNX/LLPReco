#ifndef LLPReco_DataFormats_ElectronCandidateFeatures_h
#define LLPReco_DataFormats_ElectronCandidateFeatures_h

namespace llpdnnx {
 
struct ElectronCandidateFeatures { 
    float ptrel;
    float deltaR;
    float deta;
    float dphi;
    float px;
    float py;
    float pz;
    float charge;
    float energy;
    float EtFromCaloEn;
    int isEB; 
    int isEE; 
    float ecalEnergy;
    int isPassConversionVeto;
    float convDist;
    int   convFlags;
    float convRadius;
    float hadronicOverEm;
    float ecalDrivenSeed;
    float IP2d;
    float IP2dSig;
    float IP3d;
    float IP3dSig;

    float elecSC_energy;
    float elecSC_deta;
    float elecSC_dphi;
    float elecSC_et;
    float elecSC_eSuperClusterOverP;
    float scPixCharge;
    float superClusterFbrem;

    float eSeedClusterOverP;
    float eSeedClusterOverPout;
    float eSuperClusterOverP;

    // shower shape
    float sigmaEtaEta;
    float sigmaIetaIeta;
    float sigmaIphiIphi;
    float e5x5;
    float e5x5Rel;
    float e1x5Overe5x5;
    float e2x5MaxOvere5x5;
    float r9;
    float hcalOverEcal;
    float hcalDepth1OverEcal;
    float hcalDepth2OverEcal;

    // Track-Cluster Matching Attributes
    float deltaEtaEleClusterTrackAtCalo;
    float deltaEtaSeedClusterTrackAtCalo;
    float deltaPhiSeedClusterTrackAtCalo; 
    float deltaEtaSeedClusterTrackAtVtx;
    float deltaEtaSuperClusterTrackAtVtx;
    float deltaPhiEleClusterTrackAtCalo;
    float deltaPhiSuperClusterTrackAtVtx;

    float sCseedEta;

    // electron gsf variables. 
    float EtaRel;
    float dxy;
    float dxyError;
    float dxySig;
    float dz;
    float dzError;
    float dzSig;
    int nbOfMissingHits;
    float gsfCharge;
    int   ndof;
    float chi2;
    int   numberOfBrems;
    float trackFbrem;
    float fbrem;

    // Isolation block
    float neutralHadronIso;
    float particleIso;
    float photonIso;
    float puChargedHadronIso;
    float trackIso;
    float ecalPFClusterIso;
    float hcalPFClusterIso;
    
    float pfSumPhotonEt;
    float pfSumChargedHadronPt; 
    float pfSumNeutralHadronEt;
    float pfSumPUPt;

    float dr04TkSumPt;
    float dr04EcalRecHitSumEt;
    float dr04HcalDepth1TowerSumEt;
    float dr04HcalDepth1TowerSumEtBc;
    float dr04HcalDepth2TowerSumEt;
    float dr04HcalDepth2TowerSumEtBc;
    float dr04HcalTowerSumEt;
    float dr04HcalTowerSumEtBc;
    
    
    ElectronCandidateFeatures():
        ptrel(0),
        deltaR(0),
        deta(0),
        dphi(0),
        px(0),
        py(0),
        pz(0),
        charge(0),
        energy(0),
        EtFromCaloEn(0),
        isEB(0), 
        isEE(0), 
        ecalEnergy(0),
        isPassConversionVeto(0),
        convDist(0),
        convFlags(0),
        convRadius(0),
        hadronicOverEm(0),
        ecalDrivenSeed(0),
        IP2d(0),
        IP2dSig(0),
        IP3d(0),
        IP3dSig(0),

        elecSC_energy(0),
        elecSC_deta(0),
        elecSC_dphi(0),
        elecSC_et(0),
        elecSC_eSuperClusterOverP(0),
        scPixCharge(0),
        superClusterFbrem(0),

        eSeedClusterOverP(0),
        eSeedClusterOverPout(0),
        eSuperClusterOverP(0),

        // shower shape
        sigmaEtaEta(0),
        sigmaIetaIeta(0),
        sigmaIphiIphi(0),
        e5x5(0),
        e5x5Rel(0),
        e1x5Overe5x5(0),
        e2x5MaxOvere5x5(0),
        r9(0),
        hcalOverEcal(0),
        hcalDepth1OverEcal(0),
        hcalDepth2OverEcal(0),

        // Track-Cluster Matching Attributes
        deltaEtaEleClusterTrackAtCalo(0),
        deltaEtaSeedClusterTrackAtCalo(0),
        deltaPhiSeedClusterTrackAtCalo(0), 
        deltaEtaSeedClusterTrackAtVtx(0),
        deltaEtaSuperClusterTrackAtVtx(0),
        deltaPhiEleClusterTrackAtCalo(0),
        deltaPhiSuperClusterTrackAtVtx(0),

        sCseedEta(0),

        // electron gsf variables. 
        EtaRel(0),
        dxy(0),
        dxyError(0),
        dxySig(0),
        dz(0),
        dzError(0),
        dzSig(0),
        nbOfMissingHits(0),
        gsfCharge(0),
        ndof(0),
        chi2(0),
        numberOfBrems(0),
        trackFbrem(0),
        fbrem(0),

        // Isolation block
        neutralHadronIso(0),
        particleIso(0),
        photonIso(0),
        puChargedHadronIso(0),
        trackIso(0),
        ecalPFClusterIso(0),
        hcalPFClusterIso(0),

        pfSumPhotonEt(0),
        pfSumChargedHadronPt(0), 
        pfSumNeutralHadronEt(0),
        pfSumPUPt(0),

        dr04TkSumPt(0),
        dr04EcalRecHitSumEt(0),
        dr04HcalDepth1TowerSumEt(0),
        dr04HcalDepth1TowerSumEtBc(0),
        dr04HcalDepth2TowerSumEt(0),
        dr04HcalDepth2TowerSumEtBc(0),
        dr04HcalTowerSumEt(0),
        dr04HcalTowerSumEtBc(0)
    {}
    
    bool operator<(const ElectronCandidateFeatures& other) const
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
