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

class NANOProducer : public edm::stream::EDProducer<> {
   public:
      explicit NANOProducer(const edm::ParameterSet&);
      ~NANOProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      const edm::EDGetTokenT<std::vector<reco::XTagInfo>> _tag_src;
      const edm::EDGetTokenT<std::vector<reco::LLPLabelInfo>> _label_src;
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
};

NANOProducer::NANOProducer(const edm::ParameterSet& iConfig) :
    _tag_src(consumes<std::vector<reco::XTagInfo>>(iConfig.getParameter<edm::InputTag>("srcTags"))),
    _label_src(consumes<std::vector<reco::LLPLabelInfo>>(iConfig.getParameter<edm::InputTag>("srcLabels")))
{
    produces<nanoaod::FlatTable>("global");
    produces<nanoaod::FlatTable>("csv");
    produces<nanoaod::FlatTable>("cpf");
    produces<nanoaod::FlatTable>("npf");
    produces<nanoaod::FlatTable>("sv");
    produces<nanoaod::FlatTable>("length");
    produces<nanoaod::FlatTable>("jetorigin");
}


NANOProducer::~NANOProducer()
{
}


// ------------ method called to produce the data  ------------
void
NANOProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
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

    std::vector<float> displacement;
    std::vector<float> displacement_xy;
    std::vector<float> displacement_z;
    std::vector<float> decay_angle;
    std::vector<float> betagamma;
        
    auto lengthTable = std::make_unique<nanoaod::FlatTable>(ntags, "length", false, false);
    std::vector<int> cpf_length;
    std::vector<int> npf_length;
    std::vector<int> sv_length;
    
    auto globalTable = std::make_unique<nanoaod::FlatTable>(ntags, "global", false, false);
    std::vector<float> pt;
    std::vector<float> eta;

    auto csvTable = std::make_unique<nanoaod::FlatTable>(ntags, "csv", false, false);
    std::vector<float> trackSumJetEtRatio;
    std::vector<float> trackSumJetDeltaR;
    std::vector<float> vertexCategory;
    std::vector<float> trackSip2dValAboveCharm;
    std::vector<float> trackSip2dSigAboveCharm;
    std::vector<float> trackSip3dValAboveCharm;
    std::vector<float> trackSip3dSigAboveCharm;
    std::vector<float> jetNSelectedTracks;
    std::vector<float> jetNTracksEtaRel;

    std::vector<float> cpf_trackEtaRel;
    std::vector<float> cpf_trackPtRel;
    std::vector<float> cpf_trackPPar;
    std::vector<float> cpf_trackDeltaR;
    std::vector<float> cpf_trackPParRatio;
    std::vector<float> cpf_trackPtRatio;
    std::vector<float> cpf_trackSip2dVal;
    std::vector<float> cpf_trackSip2dSig;
    std::vector<float> cpf_trackSip3dVal;
    std::vector<float> cpf_trackSip3dSig;
    std::vector<float> cpf_trackJetDistVal;
    std::vector<float> cpf_trackJetDistSig;
    std::vector<float> cpf_ptrel;
    std::vector<float> cpf_drminsv;
    std::vector<float> cpf_vertex_association;
    std::vector<float> cpf_fromPV;
    std::vector<float> cpf_puppi_weight;
    std::vector<float> cpf_track_chi2;
    std::vector<float> cpf_track_quality;
    std::vector<float> cpf_track_ndof;
    
    std::vector<float> npf_ptrel;
    std::vector<float> npf_deltaR;
    std::vector<float> npf_isGamma;
    std::vector<float> npf_hcal_fraction;
    std::vector<float> npf_drminsv;
    std::vector<float> npf_puppi_weight;
    
    std::vector<float> sv_pt;
    std::vector<float> sv_deltaR;
    std::vector<float> sv_mass;
    std::vector<float> sv_ntracks;
    std::vector<float> sv_chi2;
    std::vector<float> sv_ndof;
    std::vector<float> sv_dxy;
    std::vector<float> sv_dxysig;
    std::vector<float> sv_d3d;
    std::vector<float> sv_d3dsig;
    std::vector<float> sv_costhetasvpv;
    std::vector<float> sv_enratio;


    unsigned int ncpf_total = 0;
    unsigned int nnpf_total = 0;
    unsigned int nsv_total = 0;

    for (unsigned int itag= 0; itag < ntags; itag++) {
         const auto& features = tag_infos->at(itag).features();
         const auto& tag_info_features = features.tag_info_features;
         unsigned int ncpf = features.cpf_features.size();
         unsigned int nnpf = features.npf_features.size();
         unsigned int nsv = features.sv_features.size();
         ncpf_total += ncpf;
         nnpf_total += nnpf;
         nsv_total += nsv;
         cpf_length.push_back(ncpf);
         npf_length.push_back(nnpf);
         sv_length.push_back(nsv);

         pt.push_back(features.jet_features.pt);
         eta.push_back(features.jet_features.eta);

         trackSumJetEtRatio.push_back(tag_info_features.csv_trackSumJetEtRatio);
         trackSumJetDeltaR.push_back(tag_info_features.csv_trackSumJetDeltaR);
         trackSip2dValAboveCharm.push_back(tag_info_features.csv_trackSip2dValAboveCharm);
         trackSip2dSigAboveCharm.push_back(tag_info_features.csv_trackSip2dSigAboveCharm);
         trackSip3dValAboveCharm.push_back(tag_info_features.csv_trackSip3dValAboveCharm);
         trackSip3dSigAboveCharm.push_back(tag_info_features.csv_trackSip3dSigAboveCharm);
         jetNSelectedTracks.push_back(tag_info_features.csv_jetNSelectedTracks);
         jetNTracksEtaRel.push_back(tag_info_features.csv_jetNTracksEtaRel);
         vertexCategory.push_back(tag_info_features.csv_vertexCategory);


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
         displacement.push_back(labels.displacement);
         displacement_xy.push_back(labels.displacement_xy);
         displacement_z.push_back(labels.displacement_z);
         decay_angle.push_back(labels.decay_angle);
         betagamma.push_back(labels.betagamma);
 
    }

    auto cpfTable = std::make_unique<nanoaod::FlatTable>(ncpf_total, "cpf", false, false);
    auto npfTable = std::make_unique<nanoaod::FlatTable>(nnpf_total, "npf", false, false);
    auto svTable = std::make_unique<nanoaod::FlatTable>(nsv_total, "sv", false, false);

    for (unsigned int itag= 0; itag < ntags; itag++) {
         const auto& features = tag_infos->at(itag).features();
         auto cpf = features.cpf_features;
         auto npf = features.npf_features;
         auto sv = features.sv_features;
         unsigned int ncpf = features.cpf_features.size();
         unsigned int nnpf = features.npf_features.size();
         unsigned int nsv = features.sv_features.size();

         for (unsigned int i = 0; i < ncpf; i++){
                const auto& cpf_features = cpf.at(i);
                cpf_trackEtaRel.push_back(cpf_features.cpf_trackEtaRel);
                cpf_trackPtRel.push_back(cpf_features.cpf_trackPtRel);
                cpf_trackPPar.push_back(cpf_features.cpf_trackPPar);
                cpf_trackDeltaR.push_back(cpf_features.cpf_trackDeltaR);
                cpf_trackPParRatio.push_back(cpf_features.cpf_trackPParRatio);
                cpf_trackPtRatio.push_back(cpf_features.cpf_trackPtRatio);
                cpf_trackSip2dVal.push_back(cpf_features.cpf_trackSip2dVal);
                cpf_trackSip2dSig.push_back(cpf_features.cpf_trackSip2dSig);
                cpf_trackSip3dVal.push_back(cpf_features.cpf_trackSip3dVal);
                cpf_trackSip3dSig.push_back(cpf_features.cpf_trackSip3dSig);
                cpf_trackJetDistVal.push_back(cpf_features.cpf_trackJetDistVal);
                cpf_trackJetDistSig.push_back(cpf_features.cpf_trackJetDistSig);
                cpf_ptrel.push_back(cpf_features.cpf_ptrel);
                cpf_drminsv.push_back(cpf_features.cpf_drminsv);
                cpf_vertex_association.push_back(cpf_features.cpf_vertex_association);
                cpf_fromPV.push_back(cpf_features.cpf_fromPV);
                cpf_puppi_weight.push_back(cpf_features.cpf_puppi_weight);
                cpf_track_chi2.push_back(cpf_features.cpf_track_chi2);
                cpf_track_quality.push_back(cpf_features.cpf_track_quality);
                cpf_track_ndof.push_back(cpf_features.cpf_track_ndof);
         }
        
        for (unsigned int i = 0; i < nnpf; i++){
            const auto& npf_features = npf.at(i);
            npf_ptrel.push_back(npf_features.npf_ptrel);
            npf_deltaR.push_back(npf_features.npf_deltaR);
            npf_isGamma.push_back(npf_features.npf_isGamma);
            npf_hcal_fraction.push_back(npf_features.npf_hcal_fraction);
            npf_drminsv.push_back(npf_features.npf_drminsv);
            npf_puppi_weight.push_back(npf_features.npf_puppi_weight);
        }
        
        for (unsigned int i = 0; i < nsv; i++){
            const auto& sv_features = sv.at(i);
            sv_pt.push_back(sv_features.sv_pt);
            sv_deltaR.push_back(sv_features.sv_deltaR);
            sv_mass.push_back(sv_features.sv_mass);
            sv_ntracks.push_back(sv_features.sv_ntracks);
            sv_chi2.push_back(sv_features.sv_chi2);
            sv_ndof.push_back(sv_features.sv_ndof);
            sv_dxy.push_back(sv_features.sv_dxy);
            sv_dxysig.push_back(sv_features.sv_dxysig);
            sv_d3d.push_back(sv_features.sv_d3d);
            sv_d3dsig.push_back(sv_features.sv_d3dsig);
            sv_costhetasvpv.push_back(sv_features.sv_costhetasvpv);
            sv_enratio.push_back(sv_features.sv_enratio);
        }
    }

    globalTable->addColumn<float>("pt", pt, "global jet pt (log 10, uncorrected)", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("eta", eta, "global jet eta", nanoaod::FlatTable::FloatColumn);

    csvTable->addColumn<float>("trackSumJetEtRatio", trackSumJetEtRatio, "doc", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("trackSumJetDeltaR", trackSumJetDeltaR, "doc", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("vertexCategory", vertexCategory, "doc", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("trackSip2dValAboveCharm", trackSip2dSigAboveCharm, "doc", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("trackSip2dSigAboveCharm", trackSip2dSigAboveCharm, "doc", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("trackSip3dValAboveCharm", trackSip3dValAboveCharm, "doc", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("trackSip3dSigAboveCharm", trackSip3dSigAboveCharm, "doc", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("jetNSelectedTracks", jetNSelectedTracks, "doc", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("jetNTracksEtaRel", jetNTracksEtaRel, "doc", nanoaod::FlatTable::FloatColumn);
    
    cpfTable->addColumn<float>("trackEtaRel", cpf_trackEtaRel, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackPtRel", cpf_trackPtRel, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackPPar", cpf_trackPPar, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackDeltaR", cpf_trackDeltaR, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackPParRatio", cpf_trackPParRatio, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackPtRatio", cpf_trackPtRatio, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackSip2dVal", cpf_trackSip2dVal, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackSip2dSig", cpf_trackSip2dSig, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackSip3dVal", cpf_trackSip3dVal, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackSip3dSig", cpf_trackSip3dSig, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackJetDistVal", cpf_trackJetDistVal, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackJetDistSig", cpf_trackJetDistSig, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("ptrel", cpf_ptrel, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("drminsv", cpf_drminsv, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("vertex_association", cpf_vertex_association, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("fromPV", cpf_fromPV, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("puppi_weight", cpf_puppi_weight, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("track_chi2", cpf_track_chi2, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("track_quality", cpf_track_quality, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("track_ndof", cpf_track_ndof, "doc", nanoaod::FlatTable::FloatColumn);
    
    npfTable->addColumn<float>("ptrel", npf_ptrel, "doc", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("deltaR", npf_deltaR, "doc", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("isGamma", npf_isGamma, "doc", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("hcal_fraction", npf_hcal_fraction, "doc", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("drminsv", npf_drminsv, "doc", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("puppi_weight", npf_puppi_weight, "doc", nanoaod::FlatTable::FloatColumn);
    
    svTable->addColumn<float>("sv_pt", sv_pt, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_deltaR", sv_deltaR, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_mass", sv_mass, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_ntracks", sv_ntracks, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_chi2", sv_chi2,  "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_ndof", sv_ndof, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_dxy", sv_dxy, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_dxysig", sv_dxysig, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_d3d", sv_d3d, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_d3dsig", sv_d3dsig, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_costhetasvpv", sv_costhetasvpv, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_enratio", sv_enratio, "doc", nanoaod::FlatTable::FloatColumn);
    
    lengthTable->addColumn<int>("cpf", cpf_length, "cpf offset", nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("npf", npf_length, "npf offset", nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("sv", sv_length, "sv offset", nanoaod::FlatTable::IntColumn);

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
    
    jetOriginTable->addColumn<float>("displacement", displacement, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("displacement_xy", displacement_xy, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("displacement_z", displacement_z, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("decay_angle", decay_angle, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("betagamma", betagamma, "doc", nanoaod::FlatTable::FloatColumn);
     
    iEvent.put(std::move(globalTable), "global");
    iEvent.put(std::move(csvTable), "csv");
    iEvent.put(std::move(cpfTable), "cpf");
    iEvent.put(std::move(npfTable), "npf");
    iEvent.put(std::move(svTable), "sv");
    iEvent.put(std::move(lengthTable), "length");
    iEvent.put(std::move(jetOriginTable), "jetorigin");
}

void
NANOProducer::beginStream(edm::StreamID)
{
}

void
NANOProducer::endStream() {
}

 
void
NANOProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NANOProducer);
