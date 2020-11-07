// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "LLPReco/DataFormats/interface/XTagInfo.h"
#include "LLPReco/DataFormats/interface/LLPLabel.h"
#include "LLPReco/DataFormats/interface/LLPLabelInfo.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include <Math/Vector4D.h>




//
// class declaration
//

class NANOGenProducer : public edm::stream::EDProducer<> {
   public:
      explicit NANOGenProducer(const edm::ParameterSet&);
      ~NANOGenProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      const edm::EDGetTokenT<edm::View<pat::Jet>> _jet_src;
      const edm::EDGetTokenT<std::vector<reco::LLPLabelInfo>> _label_src;
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
};

NANOGenProducer::NANOGenProducer(const edm::ParameterSet& iConfig) :
    _jet_src(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("srcJets"))),
    _label_src(consumes<std::vector<reco::LLPLabelInfo>>(iConfig.getParameter<edm::InputTag>("srcLabels")))

{
    produces<nanoaod::FlatTable>("jetorigin");
}


NANOGenProducer::~NANOGenProducer()
{
}


// ------------ method called to produce the data  ------------
void
NANOGenProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    edm::Handle<edm::View<pat::Jet>> jets;
    iEvent.getByToken(_jet_src, jets);

    edm::Handle<std::vector<reco::LLPLabelInfo>> label_infos;
    iEvent.getByToken(_label_src, label_infos);

    std::size_t njets = jets->size();
    std::size_t ntruth = label_infos->size();


    auto jetOriginTable = std::make_unique<nanoaod::FlatTable>(ntruth, "jetorigin", false, false);



    std::vector<int> truth_jetIdx;
    std::vector<int> isPU;
    std::vector<int> isB;
    std::vector<int> isBB;
    std::vector<int> isLeptonic_B;
    std::vector<int> isLeptonic_C;
    std::vector<int> isC;
    std::vector<int> isCC;
    std::vector<int> isS;
    std::vector<int> isUD;
    std::vector<int> isG;
    std::vector<int> isPrompt_MU;
    std::vector<int> isPrompt_E;
    std::vector<int> isPrompt_TAU;
    std::vector<int> isPrompt_PHOTON;

    std::vector<int> isLLP_RAD; //no flavour match (likely from wide angle radiation)
    std::vector<int> isLLP_MU; //prompt lepton
    std::vector<int> isLLP_E; //prompt lepton
    std::vector<int> isLLP_Q; //single light quark
    std::vector<int> isLLP_QMU; //single light quark + prompt lepton
    std::vector<int> isLLP_QE; //single light quark + prompt lepton
    std::vector<int> isLLP_QQ; //double light quark
    std::vector<int> isLLP_QQMU; //double light quark + prompt lepton
    std::vector<int> isLLP_QQE; //double light quark + prompt lepton
    std::vector<int> isLLP_B; //single b/c quark
    std::vector<int> isLLP_BMU; //single b/c quark + prompt lepton
    std::vector<int> isLLP_BE; //single b/c quark + prompt lepton
    std::vector<int> isLLP_BB; //double b/c quark
    std::vector<int> isLLP_BBMU; //double b/c quark + prompt lepton
    std::vector<int> isLLP_BBE; //double b/c quark + prompt lepton
    std::vector<int> isLLP_TAU; 
    std::vector<int> isLLP_QTAU; 
    std::vector<int> isLLP_QQTAU;
    std::vector<int> isLLP_BTAU; 
    std::vector<int> isLLP_BBTAU;
    
    std::vector<int> isLLP_PHOTON;
    std::vector<int> isLLP_QPHOTON; 
    std::vector<int> isLLP_QQPHOTON;
    std::vector<int> isLLP_BPHOTON; 
    std::vector<int> isLLP_BBPHOTON;
    
    std::vector<int> isTauDecay_NO_TAU;     //no tau decay
    std::vector<int> isTauDecay_INVISIBLE;  //tau decay but not reconstructable
    std::vector<int> isTauDecay_E;          //to electron
    std::vector<int> isTauDecay_MU;         //to muon
    std::vector<int> isTauDecay_H;          //1 charged hadron
    std::vector<int> isTauDecay_H_1PI0;     //1 charged hadron + pi0(->2gamma)
    std::vector<int> isTauDecay_H_XPI0;     //1 charged hadron + 2 or more pi0(->2gamma)
    std::vector<int> isTauDecay_HHH;         //3 charged hadrons
    std::vector<int> isTauDecay_HHH_XPI0;    //3 charged hadron + 1 or more pi0(->2gamma)

    std::vector<int> isUndefined;
        
    std::vector<int> partonFlavor;
    std::vector<int> hadronFlavor;
    std::vector<int> llpId;
    std::vector<float> llp_pt;
    std::vector<float> llp_mass;

    std::vector<float> displacement;
    std::vector<float> displacement_xy;
    std::vector<float> displacement_z;
    std::vector<float> decay_angle;
    std::vector<float> betagamma;
    
    std::vector<float> matchedGenJetDeltaR;
    std::vector<float> matchedGenJetPt;
            
    for (std::size_t itag = 0; itag < ntruth; itag++) {
    	const auto& labels = label_infos->at(itag).features();
        int jetIdx = -1;
        auto base_jet_ref = label_infos->at(itag).jet();
        if (base_jet_ref.isAvailable() and base_jet_ref.isNonnull()){
        	const auto& base_jet = base_jet_ref.get();
        	for (std::size_t ijet = 0; ijet < njets; ijet++) {
            	auto jet = jets->at(ijet);
            	if (reco::deltaR(base_jet->p4(),jet.p4()) < 1e-4){
                	jetIdx = ijet;
                	break;
            	}
        	}
        } 


        truth_jetIdx.push_back(jetIdx);
        isPU.push_back(labels.type == llpdnnx::LLPLabel::Type::isPU ? 1 : 0);
        isB.push_back(labels.type == llpdnnx::LLPLabel::Type::isB ? 1 : 0);
        isBB.push_back(labels.type == llpdnnx::LLPLabel::Type::isBB ? 1 : 0);
        isLeptonic_B.push_back(labels.type == llpdnnx::LLPLabel::Type::isLeptonic_B ? 1 : 0);
        isLeptonic_C.push_back(labels.type == llpdnnx::LLPLabel::Type::isLeptonic_C ? 1 : 0);
        isC.push_back(labels.type == llpdnnx::LLPLabel::Type::isC ? 1 : 0);
        isCC.push_back(labels.type == llpdnnx::LLPLabel::Type::isCC ? 1 : 0);
        isS.push_back(labels.type == llpdnnx::LLPLabel::Type::isS ? 1 : 0);
        isUD.push_back(labels.type == llpdnnx::LLPLabel::Type::isUD ? 1 : 0);
        isG.push_back(labels.type == llpdnnx::LLPLabel::Type::isG ? 1 : 0);
        isPrompt_MU.push_back(labels.type == llpdnnx::LLPLabel::Type::isPrompt_MU ? 1 : 0);
        isPrompt_E.push_back(labels.type == llpdnnx::LLPLabel::Type::isPrompt_E ? 1 : 0);
        isPrompt_TAU.push_back(labels.type == llpdnnx::LLPLabel::Type::isPrompt_TAU ? 1 : 0);
        isPrompt_PHOTON.push_back(labels.type == llpdnnx::LLPLabel::Type::isPrompt_PHOTON ? 1 : 0);
        isLLP_RAD.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_RAD ? 1 : 0);
        isLLP_MU.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_MU ? 1 : 0);
        isLLP_E.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_E ? 1 : 0);
        isLLP_Q.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_Q ? 1 : 0);
        isLLP_QMU.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_QMU ? 1 : 0);
        isLLP_QE.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_QE ? 1 : 0);
        isLLP_QQ.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_QQ ? 1 : 0);
        isLLP_QQMU.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_QQMU ? 1 : 0);
        isLLP_QQE.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_QQE ? 1 : 0);
        isLLP_B.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_B ? 1 : 0);
        isLLP_BMU.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_BMU ? 1 : 0);
        isLLP_BE.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_BE ? 1 : 0);
        isLLP_BB.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_BB ? 1 : 0);
        isLLP_BBMU.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_BBMU ? 1 : 0);
        isLLP_BBE.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_BBE ? 1 : 0);
        isLLP_TAU.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_TAU ? 1 : 0);
        isLLP_QTAU.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_QTAU ? 1 : 0);
        isLLP_QQTAU.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_QQTAU ? 1 : 0);
        isLLP_BTAU.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_BTAU ? 1 : 0);
        isLLP_BBTAU.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_BBTAU ? 1 : 0);
        isLLP_PHOTON.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_PHOTON ? 1 : 0);
        isLLP_QPHOTON.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_QPHOTON ? 1 : 0);
        isLLP_QQPHOTON.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_QQPHOTON ? 1 : 0);
        isLLP_BPHOTON.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_BPHOTON ? 1 : 0);
        isLLP_BBPHOTON.push_back(labels.type == llpdnnx::LLPLabel::Type::isLLP_BBPHOTON ? 1 : 0);
        isUndefined.push_back(labels.type == llpdnnx::LLPLabel::Type::isUndefined ? 1 : 0);
        
        
        isTauDecay_NO_TAU.push_back(labels.tauDecay == llpdnnx::LLPLabel::TauDecay::NO_TAU ? 1 : 0);
        isTauDecay_INVISIBLE.push_back(labels.tauDecay == llpdnnx::LLPLabel::TauDecay::INVISIBLE ? 1 : 0);
        isTauDecay_E.push_back(labels.tauDecay == llpdnnx::LLPLabel::TauDecay::E ? 1 : 0);
        isTauDecay_MU.push_back(labels.tauDecay == llpdnnx::LLPLabel::TauDecay::MU ? 1 : 0);
        isTauDecay_H.push_back(labels.tauDecay == llpdnnx::LLPLabel::TauDecay::H ? 1 : 0);
        isTauDecay_H_1PI0.push_back(labels.tauDecay == llpdnnx::LLPLabel::TauDecay::H_1PI0 ? 1 : 0);
        isTauDecay_H_XPI0.push_back(labels.tauDecay == llpdnnx::LLPLabel::TauDecay::H_XPI0 ? 1 : 0);
        isTauDecay_HHH.push_back(labels.tauDecay == llpdnnx::LLPLabel::TauDecay::HHH ? 1 : 0);
        isTauDecay_HHH_XPI0.push_back(labels.tauDecay == llpdnnx::LLPLabel::TauDecay::HHH_XPI0 ? 1 : 0);


        partonFlavor.push_back(labels.partonFlavor);
        hadronFlavor.push_back(labels.hadronFlavor);
        llpId.push_back(labels.llpId);
        llp_pt.push_back(labels.llp_pt);
        llp_mass.push_back(labels.llp_mass);
        displacement.push_back(labels.displacement);
        displacement_xy.push_back(labels.displacement_xy);
        displacement_z.push_back(labels.displacement_z);
        decay_angle.push_back(labels.decay_angle);
        betagamma.push_back(labels.betagamma);
        
        matchedGenJetDeltaR.push_back(labels.matchedGenJetDeltaR);
        matchedGenJetPt.push_back(labels.matchedGenJetPt);
 
    }

    
    jetOriginTable->addColumn<int>("jetIdx", truth_jetIdx, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isPU", isPU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isB", isB, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isBB", isBB, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLeptonic_B", isLeptonic_B, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLeptonic_C", isLeptonic_C, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isC", isC, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isCC", isCC, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isS", isS, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isUD", isUD, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isG", isG, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isPrompt_MU", isPrompt_MU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isPrompt_E", isPrompt_E, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isPrompt_TAU", isPrompt_TAU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isPrompt_PHOTON", isPrompt_PHOTON, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_RAD", isLLP_RAD, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_MU", isLLP_MU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_E", isLLP_E, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_Q", isLLP_Q, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_QMU", isLLP_QMU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_QE", isLLP_QE, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_QQ", isLLP_QQ, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_QQMU", isLLP_QQMU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_QQE", isLLP_QQE, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_B", isLLP_B, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_BMU", isLLP_BMU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_BE", isLLP_BE, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_BB", isLLP_BB, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_BBMU", isLLP_BBMU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_BBE", isLLP_BBE, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_TAU", isLLP_TAU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_QTAU", isLLP_QTAU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_QQTAU", isLLP_QQTAU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_BTAU", isLLP_BTAU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_BBTAU", isLLP_BBTAU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_PHOTON", isLLP_PHOTON, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_QPHOTON", isLLP_QPHOTON, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_QQPHOTON", isLLP_QQPHOTON, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_BPHOTON", isLLP_BPHOTON, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_BBPHOTON", isLLP_BBPHOTON, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isUndefined", isUndefined, "doc", nanoaod::FlatTable::IntColumn);
    
    jetOriginTable->addColumn<int>("isTauDecay_NO_TAU", isTauDecay_NO_TAU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isTauDecay_INVISIBLE", isTauDecay_INVISIBLE, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isTauDecay_E", isTauDecay_E, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isTauDecay_MU", isTauDecay_MU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isTauDecay_H", isTauDecay_H, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isTauDecay_H_1PI0", isTauDecay_H_1PI0, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isTauDecay_H_XPI0", isTauDecay_H_XPI0, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isTauDecay_HHH", isTauDecay_HHH, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isTauDecay_HHH_XPI0", isTauDecay_HHH_XPI0, "doc", nanoaod::FlatTable::IntColumn);

    jetOriginTable->addColumn<int>("partonFlavor", partonFlavor, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("hadronFlavor", hadronFlavor, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("llpId", llpId, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<float>("llp_pt", llp_pt, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("llp_mass", llp_mass, "doc", nanoaod::FlatTable::FloatColumn);
    
    jetOriginTable->addColumn<float>("displacement", displacement, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("displacement_xy", displacement_xy, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("displacement_z", displacement_z, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("decay_angle", decay_angle, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("betagamma", betagamma, "doc", nanoaod::FlatTable::FloatColumn);
    
    jetOriginTable->addColumn<float>("matchedGenJetDeltaR", matchedGenJetDeltaR, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("matchedGenJetPt", matchedGenJetPt, "doc", nanoaod::FlatTable::FloatColumn);
     
    iEvent.put(std::move(jetOriginTable), "jetorigin");
}

void
NANOGenProducer::beginStream(edm::StreamID)
{
}

void
NANOGenProducer::endStream() {
}

 
void
NANOGenProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NANOGenProducer);

