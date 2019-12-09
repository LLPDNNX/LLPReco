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

    tensorflow::setLogging("3");
     // get threading config and build session options
    tensorflow::SessionOptions sessionOptions;
    //tensorflow::setThreading(sessionOptions, 1, "no_threads");


    // create the session using the meta graph from the cache
    _session = tensorflow::createSession(_graphDef, sessionOptions); 
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
    _input_tensors.resize(input_sizes.size());

    for (unsigned int i = 0; i < input_names_.size(); i++) {
        std::string group_name = input_names_[i];
        tensorflow::Tensor group_tensor(tensorflow::DT_FLOAT, input_sizes[i]);
        _input_tensors[i] = tensorflow::NamedTensor(group_name, group_tensor);
    }

    // Fill with zeros
    for (unsigned int i = 0; i < input_sizes.size(); i++) {
        _input_tensors[i].second.flat<float>().setZero();
    }

    // fill values of the input tensors
    for (unsigned int itag= 0; itag < ntags; itag++) {
        const auto& features = tag_infos->at(itag).features();
        //{"gen", "globalvars", "cpf", "npf", "sv"};
        // keep ctau at 1 mm for now
        auto cpf = features.cpf_features;
        unsigned int ncpf = std::min((unsigned int)cpf.size(), (unsigned int)25);
        auto npf = features.npf_features;
        unsigned int nnpf = std::min((unsigned int)npf.size(), (unsigned int)25);
        auto sv = features.sv_features;
        unsigned int nsv = std::min((unsigned int)sv.size(), (unsigned int)4);

        jet_tensor_filler(_input_tensors.at(1).second, itag, features);
        /*
        tensorflow::TTypes<float, 2>::Tensor globals = _input_tensors.at(1).second.flat_inner_dims<float>();
        for (unsigned int i = 0; i < 14; i++){
            std::cout << globals(itag, i) << " ";
        }
        */

        for (unsigned int i = 0; i < ncpf; i++){
            cpf_tensor_filler(_input_tensors.at(2).second, itag, i, cpf.at(i));
        }


        for (unsigned int i = 0; i < nnpf; i++){
            npf_tensor_filler(_input_tensors.at(3).second, itag, i, npf.at(i));
        }

        for (unsigned int i = 0; i < nsv; i++){
            sv_tensor_filler(_input_tensors.at(4).second, itag, i, sv.at(i));
        }
    }

    /*
    auto cpfs = _input_tensors.at(2).second.shaped<float, 3>({ntags, 25, 18});
    auto npfs = _input_tensors.at(3).second.shaped<float, 3>({ntags, 25, 6});
    auto svs = _input_tensors.at(4).second.shaped<float, 3>({ntags, 4, 12});
    for (unsigned int itag= 0; itag < ntags; itag++) {
        std::cout << "tag number: " << itag << std::endl;

        for (unsigned int i = 0; i < 25; i++){
            std::cout << "cand number: " << i;
            for (unsigned int j = 0; j < 18; j++){
                std::cout << cpfs(itag, i, j) << " ";
            }
            std::cout << std::endl;
        }

        for (unsigned int i = 0; i < 25; i++){
            std::cout << "cand number: " << i;
            for (unsigned int j = 0; j < 6; j++){
                std::cout << npfs(itag, i, j) << " ";
            }
            std::cout << std::endl;
        }
        for (unsigned int i = 0; i < 4; i++){
            std::cout << "cand number: " << i;
            for (unsigned int j = 0; j < 12; j++){
                std::cout << svs(itag, i, j) << " ";
            }
            std::cout << std::endl;
        }
    }
    */

    std::vector<tensorflow::Tensor> outputs;
    tensorflow::run(_session, _input_tensors, { "prediction" }, &outputs);
    tensorflow::TTypes<float, 2>::Tensor scores = outputs[0].flat_inner_dims<float>();

    for (unsigned int itag = 0; itag < ntags; itag++) {
        const auto& jet_ref = tag_infos->at(itag).jet();
        (*(output_tags))[jet_ref] = scores(itag, 4); // LLP probability
        //std::cout << scores(itag, 4) << std::endl;
    }
    iEvent.put(std::move(output_tags));
    //_session->reset();
}

void
XTagProducer::beginStream(edm::StreamID)
{
}

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
