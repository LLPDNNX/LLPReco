#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/Framework/interface/ProducerBase.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "NANOX/NANOXProducer/interface/NANOXPlugin.h"
#include "NANOX/NANOXProducer/interface/NANOXPluginFactory.h"
#include "NANOX/LLPSpecific/interface/LLPInfo.h"

#include "NANOX/DataFormats/interface/DisplacedGenVertex.h"

#include "DataFormats/Math/interface/angle.h"

#include "TRandom.h"
#include "TMath.h"

#include <iostream>

namespace nanox
{


class LLPInfoPlugin:
    public NANOXPlugin
{
    private:
        edm::EDGetTokenT<edm::View<nanox::DisplacedGenVertex>> displacedGenVertexToken_;
        std::string LLPtype;
        
        static int getHadronFlavor(const reco::Candidate& genParticle)
        {
            int absPdgId = std::abs(genParticle.pdgId());
            if (absPdgId<100)
            {
                return 0; //not a hadron
            }
            int nq3 = (absPdgId/     10)%10; //quark content
            int nq2 = (absPdgId/    100)%10; //quark content
            int nq1 = (absPdgId/   1000)%10; //quark content
            int nL  = (absPdgId/  10000)%10;
            int n   = (absPdgId/1000000)%10;
            return std::max({nq1,nq2,nq3})+n*10000+(n>0 and nL==9)*100;
        }
        
    public:
        LLPInfoPlugin(
            const std::string& name, 
            const edm::ParameterSet& pset, 
            edm::ConsumesCollector& collector,
            edm::ProducerBase& prod
        ):
            NANOXPlugin(name,pset,collector,prod),
            displacedGenVertexToken_(
                collector.consumes<edm::View<nanox::DisplacedGenVertex>>(
                    pset.getParameter<edm::InputTag>("displacedGenVertices")
                )),
            LLPtype(pset.getParameter<std::string>("LLPtype"))
        {
            prod.produces<std::vector<nanox::LLPInfo>>(name);
        }
        
        
        virtual void produce(edm::Event& event, const edm::EventSetup&) const
        {
            edm::Handle<edm::View<nanox::DisplacedGenVertex>> displacedGenVertexCollection;
            event.getByToken(displacedGenVertexToken_, displacedGenVertexCollection);
            
            std::unique_ptr<std::vector<nanox::LLPInfo>> output(new std::vector<nanox::LLPInfo>(1));
            for (size_t i = 0; i < displacedGenVertexCollection->size(); ++i)
            {
             
                nanox::LLPInfo::Data data;
                
                const nanox::DisplacedGenVertex& vertex = displacedGenVertexCollection->at(i);
                
                if (vertex.motherLongLivedParticle.isNull())
                {
                    continue;
                }
                
                auto llp = *vertex.motherLongLivedParticle.get();

                if (LLPtype == "T1qqqqLL" || LLPtype == "GMSB" ||  LLPtype == "RPV")
                {
                    if (getHadronFlavor(llp)<10000)
                    {
                        continue;
                    }
                }

                else if (LLPtype == "HToSS" || LLPtype == "HNL")
                {

                    //std::cout << "llp is : "<< llp.pdgId() <<", it coems from: " << llp.numberOfMothers() << "mothers: " << std::endl;

                    //for (unsigned int i = 0; i < llp.numberOfMothers(); i++){
                        //std::cout << i << ": " << llp.mother(i)->pdgId() << std::endl;
                    //}

                    //while (llp.numberOfMothers()>0){
                        //std::cout << " -> " << llp.mother(0)->pdgId(); 
                        //auto mother = (const reco::GenParticle*) llp.mother(0);
                        //llp = *mother;
                    //}
                    //std::cout << std::endl;
                    if (llp.numberOfMothers() == 0) continue;
                    if (LLPtype == "HToSS"){

                        if (llp.mother()->pdgId()!=25)
                        {
                            continue;
                        }
                    }

                    else if (LLPtype == "HNL"){
                        if (abs(llp.mother()->pdgId())!=24) 
                        {
                            continue;
                        }
                    }


                }


                //else if (LLPtype == "HNL")
                //{
                
                //}



                else continue;
                
                data.llp_mass = llp.mass();
                data.llp_pt = llp.pt();
                data.llp_eta = llp.eta();
                data.llp_phi = llp.phi();

                
                const reco::Candidate* lsp = nullptr;
                std::vector<int> quarks;
 
                if (LLPtype == "T1qqqqLL")
                {
               
                    for (auto genParticle: vertex.genParticles)
                    {

                        if (genParticle->numberOfDaughters()==0 and genParticle->pdgId()==1000022)
                        {
                            lsp = genParticle.get();
                        }

                       if (std::abs(genParticle->pdgId())<6 and genParticle->numberOfMothers()==1 and genParticle->mother()->pdgId()==1000021)
                        {
                            quarks.push_back(genParticle->pdgId());
                        }
                    }

                    if (quarks.size()==2 and std::abs(quarks[0])==std::abs(quarks[1]))
                    {
                        data.quarkFlavor = std::abs(quarks[0]);
                    }
 
                    if (!lsp)
                    {
                        continue;
                    }
                             
                    data.lsp_mass = lsp->mass();
                    data.lsp_pt = lsp->pt();
                    data.lsp_eta = lsp->eta();
                    data.lsp_phi = lsp->phi();

                    output->at(0).llpData.push_back(data);
                }

                else if (LLPtype == "GMSB")
                {
               
                    for (auto genParticle: vertex.genParticles)
                    {

                        if (genParticle->numberOfDaughters()==0 and genParticle->pdgId()==1000039)
                        {
                            lsp = genParticle.get();
                        }

                       if (std::abs(genParticle->pdgId())<6 and genParticle->numberOfMothers()==1 and genParticle->mother()->pdgId()==1000021)
                        {
                            quarks.push_back(genParticle->pdgId());
                        }
                    }

                    if (quarks.size()==2 and std::abs(quarks[0])==std::abs(quarks[1]))
                    {
                        data.quarkFlavor = std::abs(quarks[0]);
                    }

                    if (!lsp)
                    {
                        continue;
                    }
                             
                    data.lsp_mass = lsp->mass();
                    data.lsp_pt = lsp->pt();
                    data.lsp_eta = lsp->eta();
                    data.lsp_phi = lsp->phi();

                    output->at(0).llpData.push_back(data);
                }

                else if (LLPtype == "RPV" || LLPtype == "HToSS" || LLPtype == "HNL")
                {
               
                    output->at(0).llpData.push_back(data);
                }
                
                else continue;
            }

            event.put(std::move(output),this->name());
        }
};

}

DEFINE_EDM_PLUGIN(nanox::NANOXPluginFactory, nanox::LLPInfoPlugin, "LLPInfo");

