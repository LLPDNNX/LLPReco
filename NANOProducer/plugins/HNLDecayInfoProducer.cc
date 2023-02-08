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

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "boost/algorithm/string.hpp"

#include <vector>
#include <unordered_map>
#include <iostream>
#include <regex>


class HNLDecayInfoProducer:
    public edm::EDProducer
{
    protected:
        edm::InputTag _genParticleInputTag;
        edm::EDGetTokenT<edm::View<reco::GenParticle>> _genParticleToken;
        


    public:
        HNLDecayInfoProducer( edm::ParameterSet const & params ) :
            _genParticleInputTag(params.getParameter<edm::InputTag>("srcGenParticles")),
            _genParticleToken(consumes<edm::View<reco::GenParticle>>(_genParticleInputTag))        
        {
            produces<nanoaod::FlatTable>();
        }

        ~HNLDecayInfoProducer() override {}
        
        const reco::GenParticle* findLeptonFromTau(const reco::GenParticle* genParticle)
        {
            int absId = abs(genParticle->pdgId());
            if (absId==11 or absId==13) return genParticle;
            for (const auto daughter: genParticle->daughterRefVector())
            {
                if (daughter->isDirectPromptTauDecayProductFinalState())
                {
                    const reco::GenParticle* particle = findLeptonFromTau(daughter.get());
                    if (particle) return particle;
                }
            }
            return nullptr;
        }

        void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override 
        {
        
            edm::Handle<edm::View<reco::GenParticle>> genParticleCollection;
            iEvent.getByToken(_genParticleToken, genParticleCollection);

            const reco::GenParticle* l1 = nullptr;
            const reco::GenParticle* l2 = nullptr;

            for (const auto& genParticle: *genParticleCollection)
            {
                if (genParticle.isHardProcess())
                {
                    int absId = abs(genParticle.pdgId());
                    if (absId>=11 and absId<=16)
                    {
                        int motherId = abs(genParticle.mother(0)->pdgId());
                        if (motherId==9900012 or motherId==9990012)
                        {
                            if (l2) std::cout<<"WARNING - l2 already found"<<std::endl;
                            l2 = &genParticle;
                        }
                        else
                        {
                            if (l1) std::cout<<"WARNING - l1 already found"<<std::endl;
                            l1 = &genParticle;
                        }
                    }
                }
            }
            
            
            if (not l1)
            {
                std::cout<<"ERROR - l1 not found"<<std::endl;
            }
            if (not l2)
            {
                std::cout<<"ERROR - l2 not found"<<std::endl;
            }
            
            int l1Id = abs(l1->pdgId());
            int l2Id = abs(l2->pdgId());
            
            

            if (l1Id==15)
            {
                auto l1TauDecay = findLeptonFromTau(l1);
                if (not l1TauDecay) 
                {
                    l1Id=100; //hadronic tau
                }
                else 
                {
                    l1 = l1TauDecay;
                    l1Id=100+abs(l1->pdgId());
                }
            }
            if (l2Id==15)
            {
                auto l2TauDecay = findLeptonFromTau(l2);
                if (not l2TauDecay) 
                {
                    l2Id=100; //hadronic tau
                }
                else 
                {
                    l2 = l2TauDecay;
                    l2Id=100+abs(l2->pdgId());
                }
            }

            
            float leadingPt = 0.f;
            float trailingPt = 0.f; 
            //int lMaxId = 0;
            //int lMax = 0;
            //l1 is e or mu
            if (l1Id==11 or l1Id==13 or l1Id==111 or l1Id==113)
            {
                //lMaxId=l1Id;
                //lMax = 1;
                leadingPt = l1->pt();
                //l2 is e or mu; take pT max of l1/l2
                if (l2Id==11 or l2Id==13 or l2Id==111 or l2Id==113)
                {
                    //lMax = 2;
                    //lMaxId=l2Id;
                    
                    if (l1->pt()<l2->pt())
                    {
                        leadingPt = l2->pt();
                        trailingPt = l1->pt();
                    }
                    else
                    {
                        leadingPt = l1->pt();
                        trailingPt = l2->pt();
                    }
                }
            }
            //no l1; only l2 is e or mu
            else if (l2Id==11 or l2Id==13 or l2Id==111 or l2Id==113)
            {
                //lMax = 2;
                //lMaxId=l2Id;
                leadingPt = l2->pt();
            }
            
     
     
            auto hnlDecayTable = std::make_unique<nanoaod::FlatTable>(1,"HNLDecay",true,false);
            
            hnlDecayTable->addColumnValue<int>("l1Id",l1Id,"",nanoaod::FlatTable::IntColumn);
            hnlDecayTable->addColumnValue<float>("leadingPt",leadingPt,"",nanoaod::FlatTable::FloatColumn);
            hnlDecayTable->addColumnValue<int>("l2Id",l2Id,"",nanoaod::FlatTable::IntColumn);
            hnlDecayTable->addColumnValue<float>("trailingPt",trailingPt,"",nanoaod::FlatTable::FloatColumn);
                    
            iEvent.put(std::move(hnlDecayTable));
        }

        
        static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) 
        {
            edm::ParameterSetDescription desc;
        }

};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HNLDecayInfoProducer);

