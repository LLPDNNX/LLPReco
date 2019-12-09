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

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "FWCore/Framework/interface/makeRefToBaseProdFrom.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "LLPReco/TensorFlow/interface/tensor_fillers.h"



class XTagProducer : public edm::stream::EDProducer<> {
   public:
      explicit XTagProducer(const edm::ParameterSet&);
      ~XTagProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      const edm::EDGetTokenT<std::vector<reco::XTagInfo>> _src;
      tensorflow::GraphDef* _graphDef;
      tensorflow::Session* _session;

      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
};

XTagProducer::XTagProducer(const edm::ParameterSet& iConfig)
    : _src(consumes<std::vector<reco::XTagInfo>>(iConfig.getParameter<edm::InputTag>("src"))),
     _graphDef(tensorflow::loadGraphDef(iConfig.getParameter<edm::FileInPath>("graph_path").fullPath())),
    _session(nullptr)
{
    _session = tensorflow::createSession(_graphDef);
    produces<reco::JetTagCollection>();

}


XTagProducer::~XTagProducer()
{
      if (_session != nullptr) {
          tensorflow::closeSession(_session);
          _session = nullptr;
      }
      delete _graphDef;
      _graphDef = nullptr;
}


void
XTagProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<std::vector<reco::XTagInfo>> tag_infos;
    iEvent.getByToken(_src, tag_infos);

    // initialize output collection
    std::unique_ptr<reco::JetTagCollection> output_tags;
    if (!tag_infos->empty()) {
      auto jet_ref = tag_infos->begin()->jet();
      output_tags = std::make_unique<reco::JetTagCollection>(edm::makeRefToBaseProdFrom(jet_ref, iEvent));
    } 
    else {
      output_tags = std::make_unique<reco::JetTagCollection>();
    }

    std::vector<std::string> input_names_{"gen", "globalvars", "cpf", "npf", "sv"};
    unsigned int ntags = tag_infos->size();

    std::vector<tensorflow::TensorShape> input_sizes{
        {ntags, 1},  
        {ntags, 14},
        {ntags, 25, 18},
        {ntags, 25, 6},
        {ntags, 4, 12},
      };

    tensorflow::NamedTensorList _input_tensors;
    for (std::size_t i = 0; i < input_names_.size(); i++) {
        std::string group_name = input_names_[i];
        tensorflow::Tensor group_tensor(tensorflow::DT_FLOAT, input_sizes[i]);

        _input_tensors.push_back(tensorflow::NamedTensor(group_name, group_tensor));
    }

    // Fill with zeros
    for (std::size_t i = 0; i < input_sizes.size(); i++) {
        _input_tensors[i].second.flat<float>().setZero();
    }

    // fill values of the input tensors
    for (std::size_t itag= 0; itag < tag_infos->size(); itag++) {
        const auto& features = tag_infos->at(itag).features();
        //{"gen", "globalvars", "cpf", "npf", "sv"};
        // keep ctau at 1 mm for now
        auto cpf = features.cpf_features;
        size_t ncpf = cpf.size();
        auto npf = features.npf_features;
        size_t nnpf = npf.size();
        auto sv = features.sv_features;
        size_t nsv = sv.size();

        jet_tensor_filler(_input_tensors.at(1).second, itag, features);

        for (size_t i = 0; i < ncpf; i++){
            cpf_tensor_filler(_input_tensors.at(2).second, itag, i, cpf.at(i));
        }

        for (size_t i = 0; i < nnpf; i++){
            npf_tensor_filler(_input_tensors.at(3).second, itag, i, npf.at(i));
        }

        for (size_t i = 0; i < nsv; i++){
            sv_tensor_filler(_input_tensors.at(4).second, itag, i, sv.at(i));
        }

    }

    std::vector<tensorflow::Tensor> outputs;
    tensorflow::run(_session, _input_tensors, { "prediction" }, &outputs);

    for (unsigned int itag = 0; itag < tag_infos->size() ; itag++) {
        const auto& jet_ref = tag_infos->at(itag).jet();
        tensorflow::TTypes<float, 2>::Tensor scores = outputs[0].flat_inner_dims<float>();
        (*(output_tags))[jet_ref] = scores(itag, 0); // LLP probability?
        //std::cout << scores(itag, 0) << std::endl;
    }

    iEvent.put(std::move(output_tags));

   
    // define a tensor and fill it with range(10)
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
XTagProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
XTagProducer::endStream() {
}
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
XTagProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("pfXTagInfos"));
  desc.add<edm::FileInPath>("graph_path", edm::FileInPath("LLPReco/TensorFlow/data/da.pb"));
  descriptions.add("pfXTags", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(XTagProducer);
