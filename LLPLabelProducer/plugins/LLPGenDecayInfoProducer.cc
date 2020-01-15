// system include files
#include <memory>
#include <vector>
#include <unordered_map>
// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/Exception.h"


#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "LLPReco/DataFormats/interface/LLPGenDecayInfo.h"

#include <functional>

class LLPGenDecayInfoProducer:
    public edm::stream::EDProducer<>
    
{
    private:    
        struct LLPDecayChainCfg
        {
            std::string name;
            int llpId;
            std::vector<int> daughterIds;
        };
        
        
        static double distance(const reco::Candidate::Point& p1, const reco::Candidate::Point& p2)
        {
            return std::sqrt((p1-p2).mag2());
        }

        edm::EDGetTokenT<edm::View<reco::GenParticle>> genParticleToken_;
        
        std::vector<LLPDecayChainCfg> llpDecayChains_;
        
        virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;
            
            
         
        
    public:
        explicit LLPGenDecayInfoProducer(const edm::ParameterSet&);
        ~LLPGenDecayInfoProducer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

};



//
// constructors and destructor

//
LLPGenDecayInfoProducer::LLPGenDecayInfoProducer(const edm::ParameterSet& iConfig):
    genParticleToken_(
        consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("src"))
    )
{
    const edm::ParameterSet& decaysCfg = iConfig.getParameter<edm::ParameterSet>("decays");
    for (auto const& name: decaysCfg.getParameterNames())
    {
        const edm::ParameterSet& decayCfg = decaysCfg.getParameter<edm::ParameterSet>(name);
        
        LLPDecayChainCfg decayChainCfg;
        decayChainCfg.name = name;
        decayChainCfg.llpId = decayCfg.getParameter<int>("llpId");
        decayChainCfg.daughterIds = decayCfg.getParameter<std::vector<int>>("daughterIds");
        llpDecayChains_.push_back(decayChainCfg);
    }

    produces<std::vector<llpdnnx::LLPGenDecayInfo>>();
}


LLPGenDecayInfoProducer::~LLPGenDecayInfoProducer()
{
}


// ------------ method called to produce the data  ------------
void
LLPGenDecayInfoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<edm::View<reco::GenParticle>> genParticleCollection;
    iEvent.getByToken(genParticleToken_, genParticleCollection);
    
   
    std::unique_ptr<reco::Candidate::Point> hardInteraction(nullptr);
        
    for (unsigned int igenParticle = 0; igenParticle < genParticleCollection->size(); ++igenParticle)
    {
        const reco::GenParticle& genParticle = genParticleCollection->at(igenParticle);
        if (genParticle.isHardProcess() and genParticle.numberOfMothers()==2)
        {
            if (!hardInteraction)
            {
                hardInteraction.reset(new reco::Candidate::Point(genParticle.vertex()));
            }
            else if (distance(*hardInteraction,genParticle.vertex())>1e-10)
            {
                throw cms::Exception("PartonsFromLLPSelector: multiple hard interaction vertices found!");
            }
        }
    }
    
    
    auto output = std::make_unique<std::vector<llpdnnx::LLPGenDecayInfo>>();
    
    if (hardInteraction)
    {
        for (unsigned int igenParticle = 0; igenParticle < genParticleCollection->size(); ++igenParticle)
        {
            const reco::GenParticle& genParticle = genParticleCollection->at(igenParticle);
            
            for (auto const& decayChain: llpDecayChains_)
            {
                if (abs(genParticle.pdgId())==decayChain.llpId)
                {
                    llpdnnx::LLPGenDecayInfo llpGenDecayInfo;
                    llpGenDecayInfo.name = decayChain.name;
                    llpGenDecayInfo.llp = genParticleCollection->ptrAt(igenParticle);
            
                    for (auto const& daughter: genParticle.daughterRefVector())
                    {
                        //llp decay products need to be displaced wrt hard interaction
                        if (distance(*hardInteraction,daughter->vertex())<1e-10)
                        {
                            continue;
                        }
                        
                        if (std::find(
                            decayChain.daughterIds.cbegin(),
                            decayChain.daughterIds.cend(),
                            abs(daughter->pdgId())
                        )!=decayChain.daughterIds.cend())
                        {
                            if (daughter->pdgId()!=genParticleCollection->at(daughter.index()).pdgId())
                            {
                                throw cms::Exception("GenParticle relations not properly setup!");
                            }
                            
                            llpGenDecayInfo.decayProducts.push_back(
                                edm::Ptr<reco::GenParticle>(genParticleCollection,daughter.index())
                            );
                            
                        }
                    }
                    if (llpGenDecayInfo.decayProducts.size()>0)
                    {
                        /*
                        std::cout<<" llp decay: "<<llpGenDecayInfo.name<<", "<<llpGenDecayInfo.llp->pdgId()<<std::endl;
                        for (auto const& d: llpGenDecayInfo.decayProducts)
                        {
                            std::cout<<"   -> daughter: "<<d->pdgId()<<std::endl;
                        }
                        */
                        output->emplace_back(llpGenDecayInfo);
                    }
                }
            }
        }
    }
    
    //iEvent.put(std::move(hardInteraction));
    iEvent.put(std::move(output));
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LLPGenDecayInfoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(LLPGenDecayInfoProducer);


