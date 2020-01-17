#ifndef LLPReco_DataFormats_MuonCandidateFeatures_h
#define LLPReco_DataFormats_MuonCandidateFeatures_h


namespace llpdnnx {

class MuonCandidateFeatures {

  public:

  int  mu_jetIdx ;
  int  mu_isGlobal ;
  int  mu_isTight ;
  int  mu_isMedium ;
  int  mu_isLoose ;
  int  mu_isStandAlone ;

  float mu_pt ;
  float mu_p ;
  float mu_jetPtRel ;
  float mu_EtaRel;
  float mu_jetPtRel2 ;
  float mu_eta;
  float mu_phi;
  float mu_charge ;
  float mu_energy;
  float mu_et ;
  float mu_jetDeltaR ;
  int mu_numberOfMatchedStations ;

  float mu_2dIp ;
  float mu_2dIpSig ;
  float mu_3dIp ;
  float mu_3dIpSig ;

  float mu_absdxy ;
  float mu_absdxyError ;
  float mu_absdxySig ;
  float mu_absdz ;
  float mu_absdzError ;
  float mu_numberOfValidPixelHits;
  float mu_numberOfpixelLayersWithMeasurement ;
  float mu_numberOfstripLayersWithMeasurement ; //that does not help. needs to be discussed.

  float mu_chi2 ;
  int mu_ndof ;

  float mu_caloIso ;
  float mu_ecalIso ;
  float mu_hcalIso ;

  float mu_sumPfChHadronPt ;
  float mu_sumPfNeuHadronEt ;
  float mu_Pfpileup ;
  float mu_sumPfPhotonEt ;

  float mu_sumPfChHadronPt03 ;
  float mu_sumPfNeuHadronEt03 ;
  float mu_Pfpileup03 ;
  float mu_sumPfPhotonEt03 ;

  float mu_sumChHadronPt ;
  float mu_sumNeuHadronEt ;
  float mu_pileup ;
  float mu_sumPhotonEt ;


/*  float mu_sumChHadronPt03 ;
  float mu_sumNeuHadronEt03 ;
  float mu_pileup03 ;
  float mu_sumPhotonEt03 ; */


};


}

#endif
