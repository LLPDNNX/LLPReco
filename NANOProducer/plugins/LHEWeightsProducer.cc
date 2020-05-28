#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/NanoAOD/interface/MergeableCounterTable.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/transform.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "boost/algorithm/string.hpp"

#include <vector>
#include <unordered_map>
#include <iostream>
#include <regex>


class LHEWeightsProducer:
    public edm::EDProducer
{
    protected:
        const std::vector<edm::InputTag> lheLabel_;
        const std::vector<edm::EDGetTokenT<LHEEventProduct>> lheTag_;
        
        std::unordered_map<std::string,std::unordered_map<std::string, std::string>> weightIds_;

    public:
        LHEWeightsProducer( edm::ParameterSet const & params ) :
            lheLabel_(params.getParameter<std::vector<edm::InputTag>>("lheInfo")),
            lheTag_(edm::vector_transform(lheLabel_, [this](const edm::InputTag & tag) { return mayConsume<LHEEventProduct>(tag); }))
            
        {
            const edm::ParameterSet& pset = params.getParameter<edm::ParameterSet>("weightGroups");
            for (auto const& weightGroupName: pset.getParameterNames())
            {                
                const std::vector<std::string>& weightIdList = pset.getParameter<std::vector<std::string>>(weightGroupName);                
                
                for (size_t iweight = 0; iweight < weightIdList.size(); ++iweight)
                {
                    const auto& weightId = weightIdList[iweight];
                    weightIds_[weightGroupName][weightId] = std::to_string(iweight+1);
                }
            }
            produces<nanoaod::FlatTable>();
        }

        ~LHEWeightsProducer() override {}

        void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override 
        {
        
            edm::Handle<LHEEventProduct> lheInfo;
            for (const auto & lheTag: lheTag_)
            {
                iEvent.getByToken(lheTag, lheInfo);
                if (lheInfo.isValid()) 
                {
                    break;
                }
            }
            if (lheInfo.isValid()) 
            {
            }
            
            auto lheWeightTable = std::make_unique<nanoaod::FlatTable>(1,"LHEWeights",true,false);
            
            double normWeight = lheInfo->originalXWGTUP();
            std::vector<gen::WeightsInfo> weights = lheInfo->weights();
            for (const auto& weight: weights)
            {
                for (auto const& weightGroup: weightIds_)
                {
                    const auto weightIt = weightGroup.second.find(weight.id);
                    if (weightIt==weightGroup.second.end()) continue;

                    lheWeightTable->addColumnValue<float>(weightGroup.first+"_"+weightIt->second,weight.wgt/normWeight,"LHE weight",nanoaod::FlatTable::FloatColumn);
                    
                    //std::cout<<weight.id<<" = "<<weight.wgt<<std::endl;
                }
            }
            
            
            
            iEvent.put(std::move(lheWeightTable));
        }

        
        static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) 
        {
            edm::ParameterSetDescription desc;
        }

};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LHEWeightsProducer);

