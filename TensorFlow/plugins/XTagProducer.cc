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



class XTagProducer : public edm::stream::EDProducer<> {
   public:
      explicit XTagProducer(const edm::ParameterSet&);
      ~XTagProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      const edm::EDGetTokenT<std::vector<reco::XTagInfo>> _src;
      tensorflow::GraphDef* _graphDef;
      tensorflow::Session* _session;
      std::vector<std::pair<std::string, tensorflow::Tensor>> _inputs;

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

    std::vector<std::string> input_names_({"gen", "globalvars", "cpf", "npf", "sv"});

    std::vector<tensorflow::TensorShape> input_sizes{
        {1, 1},  
        {1, 13},
        {1, 25, 18},
        {1, 25, 6},
        {1, 4, 12},
      };

    tensorflow::NamedTensorList input_tensors;
    for (std::size_t i = 0; i < input_sizes.size(); i++) {
        input_tensors[i] =
            tensorflow::NamedTensor(input_names_[i], tensorflow::Tensor(tensorflow::DT_FLOAT, input_sizes.at(i)));
    }

    for (std::size_t i = 0; i < input_sizes.size(); i++) {
        input_tensors[i].second.flat<float>().setZero();
    }

    std::vector<tensorflow::Tensor> outputs;
    tensorflow::run(_session, input_tensors, { "output" }, &outputs);

    edm::LogInfo ("category") << " -> " << outputs[0].matrix<float>()(0, 0);

    //iEvent.put(std::move(output_tags));

   
    // define a tensor and fill it with range(10)
    /*
    tensorflow::Tensor input(tensorflow::DT_FLOAT, { 1, 10 });
    float* d = input.flat<float>().data();
    for (float i = 0; i < 10; i++, d++)
    {
        *d = i;
    }

    // define the output and run
    std::cout << "session.run" << std::endl;
    std::vector<tensorflow::Tensor> outputs;
    tensorflow::run(session_, { { "input", input } }, { "output" }, &outputs);

    // check and print the output
    std::cout << " -> " << outputs[0].matrix<float>()(0, 0) << std::endl << std::endl;


            //auto tensor_shape = _graphDef.node(inode).attr().at("shape").shape();
        //if (_graphDef.node(inode).name()==featureGroup->name())
        //{
            //foundNode = true;
            //check rank
            //if (tensor_shape.dim_size()!=(int64_t)group_shape.size())
            //{
                //throw std::runtime_error("Mismatching input rank (config: "+std::to_string(group_shape.size())+"; pb: "+std::to_string(tensor_shape.dim_size())+") for feature group '"+featureGroup->name()+"'");
            //check shape - ignore batch dim
            for (size_t idim = 1; idim < group_shape.size(); ++idim)
            {
                if ((int64_t)tensor_shape.dim(idim).size()!=group_shape[idim])
                {
                    throw std::runtime_error("Mismatching input shapes in dimension '"+std::to_string(idim+1)+"' (config: "+std::to_string(group_shape[idim])+"; pb: "+std::to_string(tensor_shape.dim(idim).size())+") for feature group '"+featureGroup->name()+"'");
                }
            }
            break;
*/
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

// ------------ method called when starting to processes a run  ------------
/*
void
XTagProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
XTagProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
XTagProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
XTagProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
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
