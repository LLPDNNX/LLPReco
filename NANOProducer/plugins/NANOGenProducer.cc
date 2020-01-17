// -*- C++ -*-
//
// Package:    LLPReco/NANOProducer
// Class:      NANOProducer
// 
/**\class NANOProducer NANOProducer.cc LLPReco/NANOProducer/plugins/NANOProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vilius Cepaitis
//         Created:  Fri, 10 Jan 2020 15:20:19 GMT
//
//


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




//
// class declaration
//

class NANOGenProducer : public edm::stream::EDProducer<> {
   public:
      explicit NANOGenProducer(const edm::ParameterSet&);
      ~NANOGenProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      const edm::EDGetTokenT<std::vector<reco::XTagInfo>> _tag_src;
      const edm::EDGetTokenT<std::vector<reco::LLPLabelInfo>> _label_src;
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
};

NANOGenProducer::NANOGenProducer(const edm::ParameterSet& iConfig) :
    _tag_src(consumes<std::vector<reco::XTagInfo>>(iConfig.getParameter<edm::InputTag>("srcTags"))),
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
    
    edm::Handle<std::vector<reco::XTagInfo>> tag_infos;
    iEvent.getByToken(_tag_src, tag_infos);

    edm::Handle<std::vector<reco::LLPLabelInfo>> label_infos;
    iEvent.getByToken(_label_src, label_infos);

    unsigned int ntags = tag_infos->size();
    unsigned int ntruth = label_infos->size();
    if (ntags != ntruth) throw cms::Exception("NANOProducer:: number of jet tags is not equal to the number of labelled jets!");
    auto jetOriginTable = std::make_unique<nanoaod::FlatTable>(ntags, "jetorigin", false, false);

    std::vector<int> label_jetIdx;
    std::vector<int> isPU;
    std::vector<int> isB;
    std::vector<int> isBB;
    std::vector<int> isGBB;
    std::vector<int> isLeptonic_B;
    std::vector<int> isLeptonic_C;
    std::vector<int> isC;
    std::vector<int> isCC;
    std::vector<int> isGCC;
    std::vector<int> isS;
    std::vector<int> isUD;
    std::vector<int> isG;

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
        
    
    for (unsigned int itag= 0; itag < ntags; itag++) {

        const auto& labels = label_infos->at(itag).features();

        unsigned int _isPU{0};
        unsigned int _isB{0};
        unsigned int _isBB{0};
        unsigned int _isGBB{0};
        unsigned int _isLeptonic_B{0};
        unsigned int _isLeptonic_C{0};
        unsigned int _isC{0};
        unsigned int _isCC{0};
        unsigned int _isGCC{0};
        unsigned int _isS{0};
        unsigned int _isUD{0};
        unsigned int _isG{0};
        unsigned int _isLLP_RAD{0};
        unsigned int _isLLP_MU{0};
        unsigned int _isLLP_E{0};
        unsigned int _isLLP_Q{0};
        unsigned int _isLLP_QMU{0};
        unsigned int _isLLP_QE{0};
        unsigned int _isLLP_QQ{0};
        unsigned int _isLLP_QQMU{0};
        unsigned int _isLLP_QQE{0};
        unsigned int _isLLP_B{0};
        unsigned int _isLLP_BMU{0};
        unsigned int _isLLP_BE{0};
        unsigned int _isLLP_BB{0};
        unsigned int _isLLP_BBMU{0};
        unsigned int _isLLP_BBE{0};
        unsigned int _isUndefined{0};

        if (labels.type == llpdnnx::LLPLabel::Type::isPU) _isPU = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isB) _isB = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isBB) _isBB = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isGBB) _isGBB = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isLeptonic_B) _isLeptonic_B = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isLeptonic_C) _isLeptonic_C = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isC) _isC = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isCC) _isCC = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isGCC) _isGCC = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isS) _isS = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isUD) _isUD = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isG) _isG = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isLLP_RAD) _isLLP_RAD = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isLLP_MU) _isLLP_MU = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isLLP_E) _isLLP_E = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isLLP_Q) _isLLP_Q = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isLLP_QMU) _isLLP_QMU = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isLLP_QE) _isLLP_QE = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isLLP_QQ) _isLLP_QQ = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isLLP_QQMU) _isLLP_QQMU = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isLLP_QQE) _isLLP_QQE = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isLLP_B) _isLLP_B = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isLLP_BMU) _isLLP_BMU = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isLLP_BE) _isLLP_BE = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isLLP_BB) _isLLP_BB = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isLLP_BBMU) _isLLP_BBMU = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isLLP_BBE) _isLLP_BBE = 1;
        if (labels.type == llpdnnx::LLPLabel::Type::isUndefined) _isUndefined = 1;


        label_jetIdx.push_back(labels.jetIdx);
        isPU.push_back(_isPU);
        isB.push_back(_isB);
        isBB.push_back(_isBB);
        isGBB.push_back(_isGBB);
        isLeptonic_B.push_back(_isLeptonic_B);
        isLeptonic_C.push_back(_isLeptonic_C);
        isC.push_back(_isC);
        isCC.push_back(_isCC);
        isGCC.push_back(_isGCC);
        isS.push_back(_isS);
        isUD.push_back(_isUD);
        isG.push_back(_isG);
        isLLP_RAD.push_back(_isLLP_RAD);
        isLLP_MU.push_back(_isLLP_MU);
        isLLP_E.push_back(_isLLP_E);
        isLLP_Q.push_back(_isLLP_Q);
        isLLP_QMU.push_back(_isLLP_QMU);
        isLLP_QE.push_back(_isLLP_QE);
        isLLP_QQ.push_back(_isLLP_QQ);
        isLLP_QQMU.push_back(_isLLP_QQMU);
        isLLP_QQE.push_back(_isLLP_QQE);
        isLLP_B.push_back(_isLLP_B);
        isLLP_BMU.push_back(_isLLP_BMU);
        isLLP_BE.push_back(_isLLP_BE);
        isLLP_BB.push_back(_isLLP_BB);
        isLLP_BBMU.push_back(_isLLP_BBMU);
        isLLP_BBE.push_back(_isLLP_BBE);
        isUndefined.push_back(_isUndefined);

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
 
    }

    
    jetOriginTable->addColumn<int>("jetIdx", label_jetIdx, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isPU", isPU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isB", isB, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isBB", isBB, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isGBB", isGBB, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLeptonic_B", isLeptonic_B, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLeptonic_C", isLeptonic_C, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isC", isC, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isCC", isCC, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isGCC", isGCC, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isS", isS, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isUD", isUD, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isG", isG, "doc", nanoaod::FlatTable::IntColumn);
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
    jetOriginTable->addColumn<int>("isUndefined", isUndefined, "doc", nanoaod::FlatTable::IntColumn);

    jetOriginTable->addColumn<int>("partonFlavor", partonFlavor, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("hadronFlavor", partonFlavor, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("llpId", partonFlavor, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<float>("llp_pt", llp_pt, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("llp_mass", llp_mass, "doc", nanoaod::FlatTable::FloatColumn);
    
    jetOriginTable->addColumn<float>("displacement", displacement, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("displacement_xy", displacement_xy, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("displacement_z", displacement_z, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("decay_angle", decay_angle, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("betagamma", betagamma, "doc", nanoaod::FlatTable::FloatColumn);
     
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

