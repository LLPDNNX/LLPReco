// system include files
#include <memory>
#include <vector>
#include <unordered_map>
// user include files

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "LLPReco/DataFormats/interface/DisplacedGenVertex.h"
#include "LLPReco/DataFormats/interface/LLPLabelInfo.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/angle.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

using llpdnnx::DisplacedGenVertex;
using llpdnnx::DisplacedGenVertexCollection;

class LLPLabelProducer:
    public edm::stream::EDProducer<>
    
{
    private:
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
    
        edm::EDGetTokenT<edm::View<pat::Jet>> jetToken_;
        edm::EDGetTokenT<edm::View<llpdnnx::DisplacedGenVertex>> displacedGenVertexToken_;

        virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;
               

    public:
        explicit LLPLabelProducer(const edm::ParameterSet&);
        ~LLPLabelProducer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

};


//
// constructors and destructor

//
LLPLabelProducer::LLPLabelProducer(const edm::ParameterSet& iConfig):
    jetToken_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("srcJets"))),
    displacedGenVertexToken_(consumes<edm::View<llpdnnx::DisplacedGenVertex>>(iConfig.getParameter<edm::InputTag>("srcVertices")))
{
    produces<reco::LLPLabelInfoCollection>();
}


LLPLabelProducer::~LLPLabelProducer()
{
}


// ------------ method called to produce the data  ------------
void
LLPLabelProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    edm::Handle<edm::View<pat::Jet>> jetCollection;
    iEvent.getByToken(jetToken_, jetCollection);
            
    edm::Handle<edm::View<llpdnnx::DisplacedGenVertex>> displacedGenVertexCollection;
    iEvent.getByToken(displacedGenVertexToken_, displacedGenVertexCollection);
    
    auto outputLLPLabelInfo = std::make_unique<reco::LLPLabelInfoCollection>();
    
    for (std::size_t ijet = 0; ijet < jetCollection->size(); ijet++) 
    {
        const pat::Jet& jet = jetCollection->at(ijet);
        edm::RefToBase<reco::Jet> jet_ref(jetCollection->refAt(ijet));
        
        llpdnnx::LLPLabel label;
        
        if (not jet.genJet())
        {
            label.type = llpdnnx::LLPLabel::Type::isPU;
        }
        else
        {
            int partonFlavor = abs(jet.partonFlavour());
            int hadronFlavor = abs(jet.hadronFlavour());
            
            label.partonFlavor = partonFlavor;
            label.hadronFlavor = hadronFlavor;
            
            unsigned int nbHadrons = jet.jetFlavourInfo().getbHadrons().size();
            unsigned int ncHadrons = jet.jetFlavourInfo().getcHadrons().size();
            
            unsigned int nbHadronsToLeptons = 0;
            unsigned int ncHadronsToLeptons = 0;
            
            unsigned int nGluons = 0; 
            
            if (jet.genParton() and (partonFlavor==5 or partonFlavor==4))
            {
                const reco::Candidate* mother = jet.genParton();
                //walk up
                while (mother->mother() and mother->mother()->pdgId()==mother->pdgId())
                {
                    mother = mother->mother();
                }
                //find gluon anchestor
                if (mother->mother() and mother->mother()->pdgId()==21)
                {
                    nGluons+=1;
                }
            }
           
            
            for (const auto* constituent: jet.genJet()->getJetConstituentsQuick())
            {   
                int absId = std::abs(constituent->pdgId());
                if (constituent->mother() and (absId==11 or absId==13))
                {
                    //account for photon/Z FSR walk up the decay tree
                    const reco::Candidate* mother = constituent->mother();
                    while (mother->mother() and mother->pdgId()==mother->mother()->pdgId())
                    {
                        mother = mother->mother();
                    }
                    int hadFlavor = getHadronFlavor(*constituent->mother());
                    if (hadFlavor==5) nbHadronsToLeptons+=1;
                    if (hadFlavor==4) ncHadronsToLeptons+=1;
                }
            }
            
            if (hadronFlavor==5)
            {
                if (nbHadronsToLeptons==0 and ncHadronsToLeptons==0)
                {
                    if (nbHadrons>1)
                    {
                        if (nGluons==0)
                        {
                            label.type = llpdnnx::LLPLabel::Type::isBB;
                        }
                        else
                        {
                            label.type = llpdnnx::LLPLabel::Type::isGBB;
                        }
                    }
                    else if (nbHadrons==1)
                    {
                        label.type = llpdnnx::LLPLabel::Type::isB;
                    }
                }
                else
                {
                    if (nbHadronsToLeptons>0)
                    {
                        label.type = llpdnnx::LLPLabel::Type::isLeptonic_B;
                    }
                    else if (nbHadronsToLeptons==0 and ncHadronsToLeptons>0)
                    {
                        label.type = llpdnnx::LLPLabel::Type::isLeptonic_C;
                    }
                }
            }
            else if (hadronFlavor==4)
            {
                if (ncHadrons>1)
                {
                    if (nGluons==0)
                    {
                        label.type = llpdnnx::LLPLabel::Type::isCC;
                    }
                    else
                    {
                        label.type = llpdnnx::LLPLabel::Type::isGCC;
                    }
                }
                else if (ncHadrons==1)
                {
                    label.type = llpdnnx::LLPLabel::Type::isC;
                }
            }
            else if (partonFlavor!=0)
            {
                if (partonFlavor==5) label.type = llpdnnx::LLPLabel::Type::isB;
                if (partonFlavor==4) label.type = llpdnnx::LLPLabel::Type::isC;
                if (partonFlavor==3) label.type = llpdnnx::LLPLabel::Type::isS;
                if (partonFlavor==2 or partonFlavor==1) label.type = llpdnnx::LLPLabel::Type::isUD;
                if (partonFlavor==21) label.type = llpdnnx::LLPLabel::Type::isG;
            }  
            else
            {
                label.type = llpdnnx::LLPLabel::Type::isUndefined;
            }
            
            
            /*
            if (displacedGenVertexCollection.product())
            {
                float dRmin = 0.4;
                //float maxDisplacement2 = 0;
                //int llpId = 0;
                //const reco::GenJet* genJetBest = nullptr;
                for (const auto& vertex: *displacedGenVertexCollection)
                {
                    for(unsigned int igenJet = 0; igenJet<vertex.genJets.size();++igenJet)
                    {
                        const reco::GenJet& genJet = vertex.genJets.at(igenJet);
                        float dRGenJets = reco::deltaR(genJet,jet);
                        //if (vertex.motherVertex.isNull()) continue;
                        //std::cout<<" - pt="<<genJet.pt()<<", dR="<<dRGenJets<<", frac="<<vertex.jetFractions[igenJet]<<std::endl;
                        //float displacement2 = (vertex.motherVertex->vertex-vertex.vertex).mag2();
                        if(dRGenJets<dRmin)
                        {
                            dRmin = dRGenJets;
                            //displacement2 = maxDisplacement2;
                            //genJetBest = &genJet;
                            if (not vertex.motherLongLivedParticle.isNull())
                            {
                                const auto &mother = *(vertex.motherLongLivedParticle);
                                if (getHadronFlavor(mother)>10000)
                                {
                                    data.fromLLP = 1;
                                }
                                //llpId = mother.pdgId();
                                
                                data.sharedVertexFraction = vertex.jetFractions[igenJet];
                                data.decay_angle = angle(genJet.p4(),mother.p4());
                                data.displacement = std::log10(std::max<float>(vertex.d3d(),1e-10));
                                data.displacement_xy = std::log10(std::max<float>(vertex.dxy(),1e-10));
                                data.displacement_z = std::log10(std::max<float>(vertex.dz(),1e-10));	
                            }
                        }
                    }		
                }
            }
            */
        }

        
        outputLLPLabelInfo->emplace_back(label,jet_ref);

    }

}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LLPLabelProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(LLPLabelProducer);

