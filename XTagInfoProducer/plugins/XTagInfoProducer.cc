/**\class XTagInfoProducer XTagInfoProducer.cc RecoLLP/XTagInfoProducer/plugins/XTagInfoProducer.cc

Description: Produces and fill in LLPDNNX features

Implementation:
    [Notes on implementation]
*/
//
// Original Author:  Vilius Cepaitis, Matthias Komm
//         Created:  Fri, 06 Dec 2019 12:09:48 GMT

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
//========= NEW ==============
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"

//============================

#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "LLPReco/DataFormats/interface/XTagInfo.h"

#include "RecoBTag/FeatureTools/interface/JetConverter.h"
#include "RecoBTag/FeatureTools/interface/ShallowTagInfoConverter.h"
#include "RecoBTag/FeatureTools/interface/SecondaryVertexConverter.h"
#include "RecoBTag/FeatureTools/interface/NeutralCandidateConverter.h"
#include "RecoBTag/FeatureTools/interface/ChargedCandidateConverter.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "RecoBTag/FeatureTools/interface/deep_helpers.h"

#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/Common/interface/Provenance.h"
#include "DataFormats/Provenance/interface/ProductProvenance.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"


#include "TVector3.h"


class XTagInfoProducer : public edm::stream::EDProducer<> {
public:
    explicit XTagInfoProducer(const edm::ParameterSet&);
    ~XTagInfoProducer();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    // ----------member data ---------------------------
    edm::EDGetTokenT<edm::View<reco::Jet>> jet_token_;
    edm::EDGetTokenT<reco::VertexCollection> vtx_token_;
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> sv_token_;
    edm::EDGetTokenT<edm::View<reco::ShallowTagInfo>> shallow_tag_info_token_;
    edm::EDGetTokenT<edm::View<reco::Candidate>> candidateToken_;
    typedef std::vector<reco::XTagInfo> XTagInfoCollection;

    edm::EDGetTokenT< pat::MuonCollection > muonsMiniAODToken_;
    edm::EDGetTokenT< pat::ElectronCollection > electronsMiniAODToken_;
 /*   edm::EDGetTokenT<edm::ValueMap<bool>>eleVetoToken_;
    edm::EDGetTokenT<edm::ValueMap<bool>>eleLooseToken_;
    edm::EDGetTokenT<edm::ValueMap<bool>>eleMediumToken_;
    edm::EDGetTokenT<edm::ValueMap<bool>>eleTightToken_;*/


};

XTagInfoProducer::XTagInfoProducer(const edm::ParameterSet& iConfig) :
jet_token_(consumes<edm::View<reco::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
vtx_token_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
sv_token_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))),
shallow_tag_info_token_(consumes<edm::View<reco::ShallowTagInfo>>(iConfig.getParameter<edm::InputTag>("shallow_tag_infos"))),
//===== NEW ======
muonsMiniAODToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
electronsMiniAODToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSrc")))

//================
{
produces<XTagInfoCollection>();
}


XTagInfoProducer::~XTagInfoProducer(){ }
void XTagInfoProducer::beginStream(edm::StreamID) { }
// ------------ method called to produce the data  ------------
void
XTagInfoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    auto output_tag_infos = std::make_unique<XTagInfoCollection>();
    edm::Handle<edm::View<reco::Jet>> jets;
    iEvent.getByToken(jet_token_, jets);

    edm::Handle<reco::VertexCollection> vtxs;
    iEvent.getByToken(vtx_token_, vtxs);

    if (vtxs->empty()) {
    // produce empty TagInfos in case no primary vertex
    iEvent.put(std::move(output_tag_infos));
    return;  // exit event
    }
    const auto& pv = vtxs->at(0);
    edm::ESHandle<TransientTrackBuilder> builder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

    edm::Handle<edm::View<reco::ShallowTagInfo>> shallow_tag_infos;
    iEvent.getByToken(shallow_tag_info_token_, shallow_tag_infos);

    edm::Handle<reco::VertexCompositePtrCandidateCollection> svs;
    iEvent.getByToken(sv_token_, svs);


//============ NEW =========== 

    edm::Handle< pat::MuonCollection > muons;
    iEvent.getByToken(muonsMiniAODToken_, muons);

    edm::Handle<pat::ElectronCollection > electrons;
    iEvent.getByToken(electronsMiniAODToken_, electrons);

//===========================

    for (std::size_t ijet = 0; ijet < jets->size(); ijet++) {
        // create data containing structure
        llpdnnx::XTagFeatures features;
        const pat::Jet& jet = jets->at(ijet); 
        edm::RefToBase<reco::Jet> jet_ref(jets, ijet);



        // Start with global jet features
        features.jet_features.pt = std::log10(jet.pt());  // uncorrected
        features.jet_features.eta = jet.eta();
        features.jet_features.mass = jet.mass();
        features.jet_features.energy = jet.energy();

        features.npv = vtxs->size();

        // Add CSV variables
        const edm::View<reco::ShallowTagInfo>& taginfos = *shallow_tag_infos;
        edm::Ptr<reco::ShallowTagInfo> match;
        // Try first by 'same index'
        if ((ijet < taginfos.size()) && (taginfos[ijet].jet() == jet_ref)) {
        match = taginfos.ptrAt(ijet);
        } else {
        // otherwise fail back to a simple search
        for (auto itTI = taginfos.begin(), edTI = taginfos.end(); itTI != edTI; ++itTI) {
            if (itTI->jet() == jet_ref) {
            match = taginfos.ptrAt(itTI - taginfos.begin());
            //break;
                }
            }
        }
        reco::ShallowTagInfo tag_info;
        if (match.isNonnull()) {
            tag_info = *match;
        }  // will be default values otherwise

        reco::TaggingVariableList vars = tag_info.taggingVariables();
        features.tag_info_features.csv_trackSumJetEtRatio = vars.get(reco::btau::trackSumJetEtRatio, -1);
        features.tag_info_features.csv_trackSumJetDeltaR = vars.get(reco::btau::trackSumJetDeltaR, -1);
        features.tag_info_features.csv_vertexCategory = vars.get(reco::btau::vertexCategory, -1);
        features.tag_info_features.csv_trackSip2dValAboveCharm = vars.get(reco::btau::trackSip2dValAboveCharm, -1);
        features.tag_info_features.csv_trackSip2dSigAboveCharm = vars.get(reco::btau::trackSip2dSigAboveCharm, -1);
        features.tag_info_features.csv_trackSip3dValAboveCharm = vars.get(reco::btau::trackSip3dValAboveCharm, -1);
        features.tag_info_features.csv_trackSip3dSigAboveCharm = vars.get(reco::btau::trackSip3dSigAboveCharm, -1);
        features.tag_info_features.csv_jetNTracksEtaRel = vars.get(reco::btau::jetNTracksEtaRel, -1);
        features.tag_info_features.csv_jetNSelectedTracks = vars.get(reco::btau::jetNSelectedTracks, -1);
                

        // fill features from secondary vertices
        for (unsigned int isv = 0; isv < svs->size(); ++isv)
        {
            const reco::VertexCompositePtrCandidate& sv = svs->at(isv);
            if (reco::deltaR(sv,jet)>0.4)
            {
                continue;
            }
            llpdnnx::SecondaryVertexFeatures sv_features;
        
            sv_features.sv_pt = std::log10(sv.pt()); 
            sv_features.sv_deltaR = reco::deltaR(sv,jet);
            sv_features.sv_mass = sv.mass();
            sv_features.sv_ntracks = sv.numberOfDaughters();
            sv_features.sv_chi2 = sv.vertexChi2();
            sv_features.sv_ndof = sv.vertexNdof();

            
            reco::Vertex::CovarianceMatrix covsv; 
            sv.fillVertexCovariance(covsv);
            reco::Vertex svtx(sv.vertex(), covsv);
        
            VertexDistanceXY distXY;
            Measurement1D distanceXY = distXY.distance(svtx, pv);
            sv_features.sv_dxy = distanceXY.value();
            sv_features.sv_dxysig = distanceXY.value()/distanceXY.error();
        
            VertexDistance3D dist3D;
            Measurement1D distance3D = dist3D.distance(svtx, pv);
            sv_features.sv_d3d = distance3D.value();
            sv_features.sv_d3dsig = distance3D.value()/distance3D.error();
            
            reco::Candidate::Vector distance(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
            sv_features.sv_costhetasvpv = sv.momentum().Unit().Dot(distance.Unit());
            sv_features.sv_enratio = sv.energy()/jet.pt();
            

            features.sv_features.emplace_back(sv_features);
        }

        std::stable_sort(features.sv_features.begin(),features.sv_features.end(),[&pv](const auto& d1, const auto& d2)
        {
            return d1.sv_dxysig>d2.sv_dxysig; //sort decreasing
        });


        // Fill track info
        for (unsigned int idaughter = 0; idaughter < jet.numberOfDaughters(); ++idaughter){
            const pat::PackedCandidate* constituent = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(idaughter));
            if ((not constituent) or constituent->charge()==0 or (not constituent->hasTrackDetails()))
                {
                    continue;
                }

            llpdnnx::ChargedCandidateFeatures cpf_features;
                
            cpf_features.cpf_ptrel = constituent->pt()/jet.pt();
                
            cpf_features.cpf_drminsv = 0.4;
            for (const auto& sv: *svs.product())
            {
                float dR = reco::deltaR(sv,*constituent);
                cpf_features.cpf_drminsv = std::min(cpf_features.cpf_drminsv,dR);
            }
                
            cpf_features.cpf_vertex_association = constituent->pvAssociationQuality();
            cpf_features.cpf_fromPV = constituent->fromPV();
            cpf_features.cpf_puppi_weight = constituent->puppiWeight();
            cpf_features.cpf_track_chi2 = constituent->pseudoTrack().chi2();
            cpf_features.cpf_track_ndof = constituent->pseudoTrack().ndof();
            cpf_features.cpf_track_quality = constituent->pseudoTrack().qualityMask();

            //if (jet.mass()<1e-10) cpf_features.cpf_jetmassdroprel = 0;
            //else cpf_features.jetmassdroprel = (jet.p4()-constituent->p4()).mass()/jet.mass();
            
            reco::TransientTrack transientTrack = builder->build(constituent->pseudoTrack());
            reco::Candidate::Vector jetDir = jet.momentum().Unit();
            GlobalVector jetRefTrackDir(jet.px(),jet.py(),jet.pz());

            Measurement1D meas_ip2d=IPTools::signedTransverseImpactParameter(transientTrack, jetRefTrackDir, pv).second;
            Measurement1D meas_ip3d=IPTools::signedImpactParameter3D(transientTrack, jetRefTrackDir, pv).second;
            Measurement1D jetdist=IPTools::jetTrackDistance(transientTrack, jetRefTrackDir, pv).second;
            reco::Candidate::Vector trackMom = constituent->pseudoTrack().momentum();
            double trackMag = std::sqrt(trackMom.Mag2());
            TVector3 trackMom3(trackMom.x(),trackMom.y(),trackMom.z());
            TVector3 jetDir3(jetDir.x(),jetDir.y(),jetDir.z());

            cpf_features.cpf_trackEtaRel=reco::btau::etaRel(jetDir, trackMom);
            cpf_features.cpf_trackPtRel=trackMom3.Perp(jetDir3);
            cpf_features.cpf_trackPPar=jetDir.Dot(trackMom);
            cpf_features.cpf_trackDeltaR=reco::deltaR(trackMom, jetDir);
            //cpf_features.cpf_trackPtRatio=cpf_features.cpf_trackPtRel / trackMag;
            cpf_features.cpf_trackPParRatio=cpf_features.cpf_trackPPar / trackMag;
            
            cpf_features.cpf_trackSip2dVal=std::abs(meas_ip2d.value());
            cpf_features.cpf_trackSip2dSig=std::abs(meas_ip2d.significance());
            cpf_features.cpf_trackSip3dVal=std::abs(meas_ip3d.value());
            cpf_features.cpf_trackSip3dSig=std::abs(meas_ip3d.significance());
            if (std::isnan(cpf_features.cpf_trackSip2dSig) || std::isnan(cpf_features.cpf_trackSip3dSig))
            {
                cpf_features.cpf_trackSip2dSig=-1.;
                cpf_features.cpf_trackSip3dSig=-1.;
            }
                
            cpf_features.cpf_trackJetDistVal = jetdist.value();
            cpf_features.cpf_trackJetDistSig = jetdist.significance();
                
            float sumPt = 0.;
            for (unsigned int jdaughter = 0; jdaughter < jet.numberOfDaughters(); ++jdaughter)
            {
                const pat::PackedCandidate* other = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(jdaughter));
                if (other and other!=constituent and reco::deltaR(*other,*constituent)<0.1)
                {
                    sumPt += other->pt();
                }
            }
            //cpf_features.cpf_relIso01 = sumPt/constituent->pt();
            
            //cpf_features.cpf_lostInnerHits = constituent->lostInnerHits(); //http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_9_4_4/doc/html/d8/d79/classpat_1_1PackedCandidate.html#ab9ef9a12f92e02fa61653ba77ee34274

            //int pdgId = std::abs(constituent->pdgId());
            //cpf_features.cpf_isLepton = pdgId==11 or pdgId==13;
        
            features.cpf_features.emplace_back(cpf_features);


// Add your code for muons and electrons here. 
//

          
              
            for(std::size_t imuon = 0 ; imuon < muons->size(); imuon++){

       //      if(jet.muonMultiplicity() == 0  ) continue ;
 
              llpdnnx::MuonCandidateFeatures mu_features ; 

              const pat::Muon& muon = muons->at(imuon);
              if(muon.isGlobalMuon() == 0) continue ; 
              if( abs(muon.eta() - constituent->eta()) < 0.01 &&  abs(muon.phi() - constituent->phi()) < 0.01 )	{
	  

                 mu_features.mu_isGlobal = muon.isGlobalMuon() ;                                   
        	 mu_features.mu_isTight = muon.isTightMuon(pv);                                     
         	 mu_features.mu_isMedium = muon.isMediumMuon();       
         	 mu_features.mu_isLoose = muon.isLooseMuon() ; 
         	 mu_features.mu_isStandAlone = muon.isStandAloneMuon() ; 

                 mu_features.mu_pt = std::log10(muon.pt());                                       
         	 mu_features.mu_p = std::log10(muon.p());  
                 mu_features.mu_jetPtRel = muon.pt()/jet.pt() ; 
         	 mu_features.mu_jetPtRel2  = muon.jetPtRatio() ; 
	         mu_features.mu_eta = muon.eta();                                                 
        	 mu_features.mu_phi = muon.phi();                                                 
	         mu_features.mu_charge = muon.charge();        
        	 mu_features.mu_energy = muon.energy();                                           
	         mu_features.mu_et = muon.et();   
                 mu_features.mu_jetDeltaR = reco::deltaR(muon ,jet) ; 
		 mu_features.mu_numberOfMatchedStations = muon.numberOfMatchedStations();

		 mu_features.mu_2dIp = muon.dB() ; 
		 mu_features.mu_2dIpSig = muon.dB()/muon.edB() ; 
		 mu_features.mu_3dIp = muon.dB(pat::Muon::PV3D) ; 
		 mu_features.mu_3dIpSig = muon.dB(pat::Muon::PV3D)/muon.edB(pat::Muon::PV3D) ;
//		 std::cout<< "1 is : "<<  mu_features.mu_2dIp << "  2 is  :  "<< mu_features.mu_2dIpSig << " 3 is :   "<< mu_features.mu_3dIp << "  4 is : "<< mu_features.mu_3dIpSig << std::endl ;

// BestTrack Info Block 
 
                 reco::Candidate::Vector muonMom = muon.bestTrack()->momentum();

                 mu_features.mu_EtaRel =reco::btau::etaRel(jetDir, muonMom);
	         mu_features.mu_absdxy = muon.bestTrack()->dxy(pv.position());
        	 mu_features.mu_absdxyError = fabs(muon.bestTrack()->dxyError()) ; 
       		 mu_features.mu_absdxySig = muon.bestTrack()->dxy(pv.position())/fabs(muon.bestTrack()->dxyError()); 
		 mu_features.mu_absdz = fabs(muon.bestTrack()->dz(pv.position())) ; 
		 mu_features.mu_absdzError = fabs(muon.bestTrack()->dzError()) ;
		 mu_features.mu_numberOfValidPixelHits = muon.bestTrack()->hitPattern().numberOfValidPixelHits();
	         mu_features.mu_numberOfpixelLayersWithMeasurement = muon.bestTrack()->hitPattern().pixelLayersWithMeasurement() ;
        	 mu_features.mu_numberOfstripLayersWithMeasurement = muon.bestTrack()->hitPattern().stripLayersWithMeasurement() ;


	         mu_features.mu_chi2 = muon.bestTrack()->normalizedChi2() ;  
		 mu_features.mu_ndof = muon.bestTrack()->ndof() ;
//		 std::cout<< "chi2 : "<< mu_features.mu_chi2 << "  mu_features.mu_ndof :  " << mu_features.mu_ndof << " mu_features.mu_numberOfValidPixelHits  " << mu_features.mu_numberOfValidPixelHits << "   mu_features.mu_numberOfpixelLayersWithMeasurement   : "<<mu_features.mu_numberOfpixelLayersWithMeasurement << "  mu_features.mu_numberOfstripLayersWithMeasurement :   "<< mu_features.mu_numberOfstripLayersWithMeasurement << std::endl ; 

// Isolation Block :
//
         	mu_features.mu_caloIso =  muon.caloIso() ; 
         	mu_features.mu_ecalIso =  muon.ecalIso() ; 
         	mu_features.mu_hcalIso =  muon.hcalIso() ;     
 
// Pf isolation : 

// Cone 0.4

	        mu_features.mu_sumPfChHadronPt  = muon.pfIsolationR04().sumChargedHadronPt;
        	mu_features.mu_sumPfNeuHadronEt  = muon.pfIsolationR04().sumNeutralHadronEt;
	        mu_features.mu_Pfpileup  = muon.pfIsolationR04().sumPUPt;
                mu_features.mu_sumPfPhotonEt = muon.pfIsolationR04().sumPhotonEt ;
          
// Cone of 0.3 

         	mu_features.mu_sumPfChHadronPt03  = muon.pfIsolationR03().sumChargedHadronPt;
         	mu_features.mu_sumPfNeuHadronEt03  = muon.pfIsolationR03().sumNeutralHadronEt;
         	mu_features.mu_Pfpileup03  = muon.pfIsolationR03().sumPUPt;
         	mu_features.mu_sumPfPhotonEt03 = muon.pfIsolationR03().sumPhotonEt ;        
         
			
		features.mu_features.emplace_back(mu_features);
  
               }

            }


    for(std::size_t ielectron = 0 ; ielectron < electrons->size(); ielectron++){
     
 //        if(jet.electronMultiplicity ()  == 0  ) continue ;

         llpdnnx::ElectronCandidateFeatures elec_features ;
     
         const pat::Electron& electron = electrons->at(ielectron);

    
         if( abs(electron.eta() - constituent->eta()) < 0.01 && abs(electron.phi() - constituent->phi()) < 0.01 ){


         elec_features.elec_pt = electron.pt() ;
	 elec_features.elec_jetPtRatio = electron.pt()/jet.pt() ; 
         elec_features.elec_p = electron.p() ; 
         elec_features.elec_eta = electron.eta() ; 
         elec_features.elec_phi = electron.phi() ; 
         elec_features.elec_charge = electron.charge() ; 
         elec_features.elec_energy = electron.energy() ; 
	 elec_features.elec_jetDeltaR = reco::deltaR(electron , jet) ; 
	 elec_features.elec_EtFromCaloEn = electron.caloEnergy() * sin(electron.p4().theta());

         elec_features.elec_isEB = electron.isEB() ;  
         elec_features.elec_isEE  = electron.isEE();
         elec_features.elec_ecalEnergy  = electron.ecalEnergy();
         elec_features.elec_isPassConversionVeto = electron.passConversionVeto();

 	 elec_features.elec_3dIP = electron.dB(pat::Electron::PV3D) ; 
	 elec_features.elec_3dIPSig = electron.dB(pat::Electron::PV3D); 
         elec_features.elec_2dIP = electron.dB() ; 
	 elec_features.elec_2dIPSig = electron.dB()/electron.edB() ; 




	 elec_features.elec_numberOfBrems  = electron.numberOfBrems () ; 
         elec_features.elec_fbrem = electron.fbrem() ; 
         elec_features.elec_e1x5 = electron.e1x5() ; 
         elec_features.elec_e2x5Max = electron.e2x5Max() ; 
         elec_features.elec_e5x5 = electron.e5x5() ;
//New : 
         elec_features.elec_hadronicOverEm = electron.hadronicOverEm() ;  
	 elec_features.elec_full5x5_sigmaIetaIeta = electron.full5x5_sigmaIetaIeta();
 	 elec_features.elec_deltaEtaEleClusterTrackAtCalo  = electron.deltaEtaEleClusterTrackAtCalo();
	 elec_features.elec_deltaEtaSeedClusterTrackAtCalo = electron.deltaEtaSeedClusterTrackAtCalo () ; 
	 elec_features.elec_deltaEtaSeedClusterTrackAtVtx = electron.deltaEtaSeedClusterTrackAtVtx();

	 elec_features.elec_deltaEtaSuperClusterTrackAtVtx = electron.deltaEtaSuperClusterTrackAtVtx() ;  
	 elec_features.elec_deltaPhiEleClusterTrackAtCalo = electron.deltaPhiEleClusterTrackAtCalo() ; 
	 elec_features.elec_deltaPhiSeedClusterTrackAtCalo = electron.deltaPhiSeedClusterTrackAtCalo() ; 
	 elec_features.elec_deltaPhiSuperClusterTrackAtVtx = electron.deltaPhiSuperClusterTrackAtVtx ()  ;

//GSF -> Gaussian Sum Filter
 
	 elec_features.elec_dxy = electron.gsfTrack()->dxy(pv.position()) ; 
	 elec_features.elec_dz = electron.gsfTrack()->dz(pv.position()) ;
	 elec_features.elec_nbOfMissingHits = electron.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) ; 
         elec_features.elec_gsfCharge = electron.gsfTrack()->charge() ;
 
// Super Cluster variables.
//
  	 elec_features.elecSC_energy = electron.superCluster()->energy() ; 
	 elec_features.elecSC_eta = electron.superCluster()->eta();
 	 elec_features.elecSC_phi = electron.superCluster()->phi();
         elec_features.elecSC_et = electron.superCluster()->energy() * sin(electron.p4().theta());
         elec_features.elecSC_eSuperClusterOverP  = electron.eSuperClusterOverP();

	
	 elec_features.elec_scE1x5 = electron.scE1x5 () ; 
	 elec_features.elec_scE2x5Max  = electron.scE2x5Max() ; 
	 elec_features.elec_scE5x5 = electron.scE5x5 () ; 
	 elec_features.elec_scPixCharge = electron.scPixCharge() ; 
	 elec_features.elec_scSigmaEtaEta = electron.scSigmaEtaEta() ; 
	 elec_features.elec_scSigmaIEtaIEta = electron.scSigmaIEtaIEta() ;
	 elec_features.elec_superClusterFbrem = electron.superClusterFbrem() ; 

// electron Isolation 
//
	elec_features.elec_hcalPFClusterIso = electron.hcalPFClusterIso() ; 
        elec_features.elec_dr03TkSumPt = electron.dr03TkSumPt() ;
 
        elec_features.elec_hcalDepth1OverEcal = electron.hcalDepth1OverEcal() ; 

	elec_features.elec_hcalDepth2OverEcal = electron.hcalDepth2OverEcal() ; 

	elec_features.elec_dr03HcalDepth2TowerSumEt = electron.dr03HcalDepth2TowerSumEt() ; 

	elec_features.elec_hcalDepth2TowerSumEtNoVeto = electron.isolationVariables03().hcalDepth2TowerSumEt ; 

        elec_features.elec_hcalDepth1TowerSumEtNoVeto = electron.isolationVariables03().hcalDepth1TowerSumEt ; 

        elec_features.elec_dr03EcalRecHitSumEt = electron.dr03EcalRecHitSumEt() ; 

	elec_features.elec_dr03HcalDepth1TowerSumEt = electron.dr03HcalDepth1TowerSumEt() ; 

	elec_features.elec_dr03HcalDepth1TowerSumEtBc = electron.dr03HcalDepth1TowerSumEtBc(); 

        elec_features.elec_pfSumPhotonEt = electron.pfIsolationVariables().sumPhotonEt ; 

	elec_features.elec_pfSumChargedHadronPt = electron.pfIsolationVariables().sumChargedHadronPt; 

        elec_features.elec_pfSumNeutralHadronEt = electron.pfIsolationVariables().sumNeutralHadronEt ; 

        elec_features.elec_pfSumPUPt = electron.pfIsolationVariables().sumPUPt ; 	
     
	features.elec_features.emplace_back(elec_features);
        }

   }

 
        }
        std::stable_sort(features.cpf_features.begin(),features.cpf_features.end(),[](const auto& d1, const auto& d2)
            {
                if (!std::isnan(d1.cpf_trackSip2dSig) and !std::isinf(d1.cpf_trackSip2dSig))
                {
                    if (!std::isnan(d2.cpf_trackSip2dSig) and !std::isinf(d2.cpf_trackSip2dSig))
                    {
                        if (std::fabs(d1.cpf_drminsv-d2.cpf_drminsv)>std::numeric_limits<float>::epsilon())
                        {
                            return std::fabs(d1.cpf_trackSip2dSig)>std::fabs(d2.cpf_trackSip2dSig); //sort decreasing
                        }
                    }
                    else
                    {
                        return true;
                    }
                }
                else if (!std::isnan(d2.cpf_trackSip2dSig) and !std::isinf(d2.cpf_trackSip2dSig))
                {
                    return false;
                }
                if (!std::isnan(d1.cpf_drminsv) and !std::isinf(d1.cpf_drminsv))
                {
                    if (!std::isnan(d2.cpf_drminsv) and !std::isinf(d2.cpf_drminsv))
                    {
                        if (std::fabs(d1.cpf_drminsv-d2.cpf_drminsv)>std::numeric_limits<float>::epsilon())
                        {
                            return d1.cpf_drminsv<d2.cpf_drminsv; //sort increasing
                        }
                    }
                    else
                    {
                        return false;
                    }
                }
                else if (!std::isnan(d2.cpf_drminsv) and !std::isinf(d2.cpf_drminsv))
                {
                    return true;
                }
                
                if (!std::isnan(d1.cpf_ptrel) and !std::isinf(d1.cpf_ptrel))
                {
                    if (!std::isnan(d2.cpf_ptrel) and !std::isinf(d2.cpf_ptrel))
                    {
                        return d1.cpf_ptrel>d2.cpf_ptrel; //sort decreasing
                    }
                    else
                    {
                        return true;
                    }
                }
                else if (!std::isnan(d2.cpf_ptrel) and !std::isinf(d2.cpf_ptrel))
                {
                    return false;
                }
                return false;
            });
        // Fill neutral hadron info
        for (unsigned int idaughter = 0; idaughter < jet.numberOfDaughters(); ++idaughter){
            const pat::PackedCandidate* constituent = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(idaughter));
            if ((not constituent) or constituent->charge()!=0)
            {
                continue;
            }
            llpdnnx::NeutralCandidateFeatures npf_features;
        
            npf_features.npf_ptrel = constituent->pt()/jet.pt();
            npf_features.npf_puppi_weight = constituent->puppiWeight();
            npf_features.npf_deltaR = reco::deltaR(*constituent,jet);
            npf_features.npf_isGamma = fabs(constituent->pdgId())==22;
            npf_features.npf_hcal_fraction = constituent->hcalFraction();
            
            npf_features.npf_drminsv = 0.4;
            for (const auto& sv: *svs.product())
            {
                float dR = reco::deltaR(sv,*constituent);
                npf_features.npf_drminsv = std::min(npf_features.npf_drminsv,dR);
            }
            
            //if (jet.mass()<1e-10) npf_features.npf._etmassdroprel = -9;
            //else npf_features.npf_jetmassdroprel = (jet.p4()- constituent->p4()).mass()/jet.mass();

            /*
            float sumPt = 0.;
            for (unsigned int jdaughter = 0; jdaughter < jet.numberOfDaughters(); ++jdaughter)
            {
                const pat::PackedCandidate* other = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(jdaughter));
                if (other and other!=constituent and reco::deltaR(*other,*constituent)<0.1)
                {
                    sumPt += other->pt();
                }
            }
            npf_features.npf_relIso01 = sumPt/constituent->pt();
            */
            features.npf_features.emplace_back(npf_features);
        }
        std::stable_sort(features.npf_features.begin(),features.npf_features.end(),[](const auto& d1, const auto& d2)
            {
                if (!std::isnan(d1.npf_drminsv) and !std::isinf(d1.npf_drminsv))
                {
                    if (!std::isnan(d2.npf_drminsv) and !std::isinf(d2.npf_drminsv))
                    {
                        if (std::fabs(d1.npf_drminsv-d2.npf_drminsv)>std::numeric_limits<float>::epsilon())
                        {
                            return d1.npf_drminsv<d2.npf_drminsv; //sort increasing
                        }
                    }
                    else
                    {
                        return false;
                    }
                }
                else if (!std::isnan(d2.npf_drminsv) and !std::isinf(d2.npf_drminsv))
                {
                    return true;
                }
                    
                if (!std::isnan(d1.npf_ptrel) and !std::isinf(d1.npf_ptrel))
                    {
                    if (!std::isnan(d2.npf_ptrel) and !std::isinf(d2.npf_ptrel))
                    {
                        return d1.npf_ptrel>d2.npf_ptrel; //sort decreasing
                    }
                    else
                    {
                        return true;
                    }
                }
                else if (!std::isnan(d2.npf_ptrel) and !std::isinf(d2.npf_ptrel))
                {
                    return false;
                }
                return false;
            });
    output_tag_infos->emplace_back(features, jet_ref);
  }

  iEvent.put(std::move(output_tag_infos));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void XTagInfoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
edm::ParameterSetDescription desc;
desc.add<edm::InputTag>("jets", edm::InputTag("ak4PFJetsCHS"));
desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("inclusiveCandidateSecondaryVertices"));
desc.add<edm::InputTag>("shallow_tag_infos", edm::InputTag("pfDeepCSVTagInfos"));
}
void XTagInfoProducer::endStream() {};

//define this as a plug-in
DEFINE_FWK_MODULE(XTagInfoProducer);
