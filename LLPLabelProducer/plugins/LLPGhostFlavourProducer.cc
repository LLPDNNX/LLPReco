// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "PhysicsTools/JetMCUtils/interface/CandMCTag.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "LLPReco/DataFormats/interface/LLPGenDecayInfo.h"
#include "LLPReco/DataFormats/interface/LLPGhostFlavourInfo.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"


//this class will be attached to the pseudojet to identify the ghost in the clustered jets
class GhostInfo: 
    public fastjet::PseudoJet::UserInfoBase
{
    public:
        GhostInfo(
            size_t decayIndex,
            size_t decayProductIndex
        ):
            decayIndex_(decayIndex),
            decayProductIndex_(decayProductIndex)
        {
        }
        
        inline size_t decayIndex() const
        { 
            return decayIndex_;
        }
        
        inline size_t decayProductIndex() const
        { 
            return decayProductIndex_;
        }

    protected:
        size_t decayIndex_;
        size_t decayProductIndex_;
};

class LLPGhostFlavourProducer: 
    public edm::EDProducer
{
    public:
        explicit LLPGhostFlavourProducer(const edm::ParameterSet&);
        ~LLPGhostFlavourProducer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void produce(edm::Event&, const edm::EventSetup&);

        void insertGhosts(const edm::Handle<reco::GenParticleRefVector>& particles,
                    const double ghostRescaling,
                    std::vector<fastjet::PseudoJet>& constituents);

        void matchReclusteredJets(const edm::Handle<edm::View<pat::Jet> >& jets,
                            const std::vector<fastjet::PseudoJet>& matchedJets,
                            std::vector<int>& matchedIndices);


        const edm::EDGetTokenT<edm::View<pat::Jet>>  jetsToken_;   
        const edm::EDGetTokenT<edm::View<llpdnnx::LLPGenDecayInfo>> llpGenDecayInfoToken_;

        const std::string jetAlgorithm_;
        const double rParam_;
        const double ghostRescaling_;
        const double relPtTolerance_;

        fastjet::JetDefinition fjJetDefinition_;
};


template<typename HANDLE, typename TYPE>
void writeValueMap(edm::Event &out, const HANDLE& handle, const std::vector<TYPE> values, const std::string &name)  
{
     typedef edm::ValueMap<TYPE> Map;
     std::unique_ptr<Map> map(new Map());
     typename Map::Filler filler(*map);
     filler.insert(handle, values.begin(), values.end());
     filler.fill();
     out.put(std::move(map), name);
}


void getAllDaughers(const reco::GenParticle* p, std::vector<const reco::GenParticle*>& daughters)
{
    for (const auto& d: p->daughterRefVector ())
    {
        const reco::GenParticle* dcast = dynamic_cast<const reco::GenParticle*>(d.get());
        if (not d)
        {
            std::cout<<"ERROR while casting genparticle"<<std::endl;
        }
        if (dcast->	numberOfDaughters ()==0)
        {
            daughters.push_back(dcast);
        }
        else
        {
            getAllDaughers(dcast, daughters);
        }
    }
}


LLPGhostFlavourProducer::LLPGhostFlavourProducer(const edm::ParameterSet& iConfig) :
    jetsToken_(consumes<edm::View<pat::Jet> >( iConfig.getParameter<edm::InputTag>("srcJets")) ),
    llpGenDecayInfoToken_(consumes<edm::View<llpdnnx::LLPGenDecayInfo>>( iConfig.getParameter<edm::InputTag>("srcDecayInfo") )),
    jetAlgorithm_(iConfig.getParameter<std::string>("jetAlgorithm")),
    rParam_(iConfig.getParameter<double>("rParam")),
    ghostRescaling_(iConfig.exists("ghostRescaling") ? iConfig.getParameter<double>("ghostRescaling") : 1e-18),
    relPtTolerance_(iConfig.exists("relPtTolerance") ? iConfig.getParameter<double>("relPtTolerance") : 1e-04) // 0.1% relative difference in Pt should be sufficient to detect possible misconfigurations

{
    produces<edm::ValueMap<llpdnnx::LLPGhostFlavourInfo>>();

    // set jet algorithm
    if (jetAlgorithm_=="Kt")
    {
        fjJetDefinition_= fastjet::JetDefinition(fastjet::kt_algorithm, rParam_);
    }
    else if (jetAlgorithm_=="CambridgeAachen")
    {
        fjJetDefinition_= fastjet::JetDefinition(fastjet::cambridge_algorithm, rParam_);
    }
    else if (jetAlgorithm_=="AntiKt")
    {
        fjJetDefinition_= fastjet::JetDefinition(fastjet::antikt_algorithm, rParam_);
    }
    else
    {
        throw cms::Exception("InvalidJetAlgorithm") << "Jet clustering algorithm is invalid: " << jetAlgorithm_ << ", use CambridgeAachen | Kt | AntiKt" << std::endl;
    }
}


LLPGhostFlavourProducer::~LLPGhostFlavourProducer()
{
}


void
LLPGhostFlavourProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<edm::View<pat::Jet> > jetCollection;
    iEvent.getByToken(jetsToken_, jetCollection);

    edm::Handle<edm::View<llpdnnx::LLPGenDecayInfo>> llpGenDecayInfoCollection;
    iEvent.getByToken(llpGenDecayInfoToken_, llpGenDecayInfoCollection);

    // vector of constituents for reclustering jets and "ghosts"
    std::vector<fastjet::PseudoJet> fjInputs;
    // loop over all input jets and collect all their constituents
    for (auto const& jet: *jetCollection)
    {
        for(unsigned int iconstituent = 0; iconstituent < jet.numberOfDaughters(); ++iconstituent)
        {
            const reco::Candidate* constituent = jet.daughter(iconstituent);
            if(constituent->pt() < 1e-10)
            {
                edm::LogWarning("NullTransverseMomentum") << "dropping input candidate with pt<1e-10";
                continue;
            }
            fjInputs.push_back(fastjet::PseudoJet(constituent->px(),constituent->py(),constituent->pz(),constituent->energy()));
        }
    }
    
    //insertGhosts
    for (size_t idecay = 0; idecay < llpGenDecayInfoCollection->size(); ++idecay)
    {
        auto const& decayInfo = llpGenDecayInfoCollection->at(idecay);
        for (size_t idaughter = 0; idaughter < decayInfo.decayProducts.size(); ++idaughter)
        {
            auto const& daughter = decayInfo.decayProducts[idaughter];
            fastjet::PseudoJet pseudojet(daughter->px(),daughter->py(),daughter->pz(),daughter->energy());
            pseudojet*=ghostRescaling_; // rescale particle momentum
            pseudojet.set_user_info(new GhostInfo(idecay,idaughter));
            fjInputs.push_back(pseudojet);
        }
    }
    
   
    // define jet clustering sequence
    fastjet::ClusterSequence fjClusterSequence(fjInputs, fjJetDefinition_);
    std::vector<fastjet::PseudoJet> reclusteredJets = fastjet::sorted_by_pt(
        fjClusterSequence.inclusive_jets(ghostRescaling_*0.001)
    );

    if (reclusteredJets.size() < jetCollection->size())
    {
        edm::LogError("TooFewReclusteredJets") 
            << "There are fewer reclustered (" 
            << reclusteredJets.size() 
            << ") than original jets (" 
            << jetCollection->size() 
            << "). Please check that the jet algorithm and jet size match those used for the original jet collection.";
    }
    // match reclustered and original jets
    std::vector<int> reclusteredIndices;
    matchReclusteredJets(jetCollection,reclusteredJets,reclusteredIndices);

    // determine jet flavour
    std::vector<llpdnnx::LLPGhostFlavourInfo> ghostFlavourInfo(jetCollection->size());
    for(size_t ijet = 0; ijet < jetCollection->size(); ++ijet)
    {
        if (reclusteredIndices[ijet]<0) //jet is not matched
        {
            continue;
        }
        
        const pat::Jet& recoJet = jetCollection->at(ijet);
        const fastjet::PseudoJet& reclusteredJet = reclusteredJets[reclusteredIndices[ijet]];
        
        if(recoJet.pt() < 1e-10 )
        {
            edm::LogWarning("NullTransverseMomentum") 
                << "The original jet " << ijet
                << " has Pt=0. This is not expected so the jet will be skipped.";
        }
        else
        {

            // since the "ghosts" are extremely soft, the configuration and ordering of the reclustered and original jets should in principle stay the same
            double ptRel = std::fabs(reclusteredJet.pt() - recoJet.correctedP4("Uncorrected").pt()) / recoJet.correctedP4("Uncorrected").pt();
            if (ptRel > relPtTolerance_)
            {
                edm::LogError("JetPtMismatch") 
                    << "The reclustered and original jet " << ijet << " have different Pt's (" 
                    << reclusteredJet.pt() << " vs " 
                    << recoJet.correctedP4("Uncorrected").pt() << " GeV, respectively)";
            }

            // loop over jet constituents and try to find "ghosts"
            for(auto const& constituent: reclusteredJet.constituents())
            {
                if (not constituent.has_user_info()) continue; // skip if not a "ghost"
                GhostInfo ghostInfo = constituent.user_info<GhostInfo>();
                
                
                ghostFlavourInfo[ijet].llpFlavours.emplace_back(
                    llpGenDecayInfoCollection->ptrAt(ghostInfo.decayIndex()),
                    llpGenDecayInfoCollection->at(ghostInfo.decayIndex()).decayProducts[ghostInfo.decayProductIndex()]
                );
            }
        }
    }
    
    writeValueMap(iEvent,jetCollection,ghostFlavourInfo,"");
}


void
LLPGhostFlavourProducer::matchReclusteredJets(
    const edm::Handle<edm::View<pat::Jet> >& jetCollection,
    const std::vector<fastjet::PseudoJet>& reclusteredJets,
    std::vector<int>& matchedIndices
)
{
    std::vector<bool> matchedLocks(reclusteredJets.size(),false);

    for (auto const& jet: *jetCollection)
    {
        double matchedDR2 = 1e9;
        int matchedIdx = -1;

        for(size_t rj=0; rj<reclusteredJets.size(); ++rj)
        {
            if( matchedLocks.at(rj) ) continue; // skip jets that have already been matched

            double tempDR2 = reco::deltaR2(
                jet.rapidity(), 
                jet.phi(), 
                reclusteredJets.at(rj).rapidity(), 
                reclusteredJets.at(rj).phi_std()
            );
            if( tempDR2 < matchedDR2 )
            {
                matchedDR2 = tempDR2;
                matchedIdx = rj;
            }
        }

        if( matchedIdx>=0 )
        {
            if ( matchedDR2 > rParam_*rParam_ )
            {
                edm::LogError("JetMatchingFailed") 
                    << "Matched reclustered jet and original jet are separated by dR=" 
                    << sqrt(matchedDR2) 
                    << " which is greater than the jet size R=" 
                    << rParam_ << ".\n"
                    << "This is not expected so please check that the jet algorithm and jet size match those used for the original jet collection.";
            }
            else
            {
                matchedLocks.at(matchedIdx) = true;
            }
        }
        else
        {
            edm::LogError("JetMatchingFailed") 
                << "Matching reclustered to original jets failed. "
                << "Please check that the jet algorithm and jet size match those used for the original jet collection.";
        }
        matchedIndices.push_back(matchedIdx);
    }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LLPGhostFlavourProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LLPGhostFlavourProducer);

