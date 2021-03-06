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

#include "LLPReco/DataFormats/interface/LLPGhostFlavourInfo.h"

#include "TH1.h"


using llpdnnx::DisplacedGenVertex;
using llpdnnx::DisplacedGenVertexCollection;

class LLPLabelOldProducer:
    public edm::stream::EDProducer<>
    
{
    private:
        struct LLPDecay
        {
            const reco::GenParticle* llp;
            std::vector<const reco::GenParticle*> decayProducts;
            
            LLPDecay():
                llp(nullptr)
            {
            }
        };

        //TH1D *hist;
        //edm::Service<TFileService> fs;
    
    
        static int getHadronFlavor(const reco::Candidate& genParticle)
        {
            int absPdgId = std::abs(genParticle.pdgId());
            if (absPdgId<100)
            {
                if (absPdgId<6) return absPdgId; //parton
                return 0; //not a hadron
            }
            int nq3 = (absPdgId/     10)%10; //quark content
            int nq2 = (absPdgId/    100)%10; //quark content
            int nq1 = (absPdgId/   1000)%10; //quark content
            int nL  = (absPdgId/  10000)%10;
            int n   = (absPdgId/1000000)%10;
            return std::max({nq1,nq2,nq3})+n*10000+(n>0 and nL==9)*100;
        }
        double tauPtThreshold_;
        double quarkPtThreshold_;
        double bPtThreshold_;
        double muonPtThreshold_;
        double electronPtThreshold_;
        
        std::unordered_map<llpdnnx::LLPLabel::Type,int> jetsPerLabel_;
    
        edm::EDGetTokenT<edm::View<pat::Jet>> jetToken_;
        edm::EDGetTokenT<edm::View<llpdnnx::DisplacedGenVertex>> displacedGenVertexToken_;
        edm::EDGetTokenT<edm::ValueMap<llpdnnx::LLPGhostFlavourInfo>> llpFlavourInfoToken_;


        virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

    public:
        explicit LLPLabelOldProducer(const edm::ParameterSet&);
        ~LLPLabelOldProducer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

};

LLPLabelOldProducer::LLPLabelOldProducer(const edm::ParameterSet& iConfig):
    tauPtThreshold_(iConfig.getParameter<double>("tauPtThreshold")),
    quarkPtThreshold_(iConfig.getParameter<double>("quarkPtThreshold")),
    bPtThreshold_(iConfig.getParameter<double>("bPtThreshold")),
    muonPtThreshold_(iConfig.getParameter<double>("muonPtThreshold")),
    electronPtThreshold_(iConfig.getParameter<double>("electronPtThreshold")),
    jetToken_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("srcJets"))),
    displacedGenVertexToken_(consumes<edm::View<llpdnnx::DisplacedGenVertex>>(iConfig.getParameter<edm::InputTag>("srcVertices"))),
    llpFlavourInfoToken_(consumes<edm::ValueMap<llpdnnx::LLPGhostFlavourInfo>>(iConfig.getParameter<edm::InputTag>("srcFlavourInfo")))
{
    produces<reco::LLPLabelInfoCollection>();
    //hist = fs->make<TH1D>("ptfrac" , "ptfrac" , 100 , 0 , 1. );
}


LLPLabelOldProducer::~LLPLabelOldProducer()
{
    for (const auto& labelCountPair: jetsPerLabel_)
    {
        std::cout<<llpdnnx::LLPLabel::typeToString(labelCountPair.first)<<": "<<labelCountPair.second<<std::endl;
    }
}



// ------------ method called to produce the data  ------------
void
LLPLabelOldProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    edm::Handle<edm::View<pat::Jet>> jetCollection;
    iEvent.getByToken(jetToken_, jetCollection);
          
    edm::Handle<edm::View<llpdnnx::DisplacedGenVertex>> displacedGenVertexCollection;
    iEvent.getByToken(displacedGenVertexToken_, displacedGenVertexCollection);
    
    edm::Handle<edm::ValueMap<llpdnnx::LLPGhostFlavourInfo>> llpGhostFlavourInfoMap;
    iEvent.getByToken(llpFlavourInfoToken_, llpGhostFlavourInfoMap);
    
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
            label.matchedGenJetDeltaR = reco::deltaR(jet.p4(),jet.genJet()->p4());
            label.matchedGenJetPt = jet.genJet()->pt();
            
            unsigned int nbHadrons = jet.jetFlavourInfo().getbHadrons().size();
            unsigned int ncHadrons = jet.jetFlavourInfo().getcHadrons().size();

            unsigned int nbHadronsToLeptons = 0;
            unsigned int ncHadronsToLeptons = 0;
            
            unsigned int nPromptElectrons = 0;
            unsigned int nPromptMuons = 0;

            unsigned int nGluons = 0; 

            float tauPtFrac = 0;

            float promptPtThreshold = 0.6;
            
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

            for (unsigned int iConst = 0; iConst < jet.genJet()->numberOfDaughters(); iConst++)
            {
                const reco::Candidate* constituent = jet.genJet()->daughter(iConst);
                const pat::PackedGenParticle* packedConstituent = dynamic_cast<const pat::PackedGenParticle*>(constituent);
                unsigned int absId = std::abs(constituent->pdgId());
                if (not (constituent->mother()))
                {
                    continue;
                }

                if (absId == 11 or absId == 13)
                {
                    if (packedConstituent->isPromptFinalState() or packedConstituent->isDirectPromptTauDecayProductFinalState())
                    {
                        int hadFlavor = getHadronFlavor(*constituent->mother());
                        if (hadFlavor==5) nbHadronsToLeptons+=1;
                        if (hadFlavor==4) ncHadronsToLeptons+=1;                   

                        float ptFrac = constituent->pt()/jet.genJet()->pt();
                        if (ptFrac > promptPtThreshold)
                        {
                            if (absId == 13){
                                nPromptMuons += 1;
                            }
                            else if (absId == 11){
                                nPromptElectrons += 1;
                            }           
                        }
                    }
                }

                else if (packedConstituent->isDirectPromptTauDecayProductFinalState())
                {
                    tauPtFrac += constituent->pt()/jet.genJet()->pt();
                }
            }


            if (tauPtFrac > promptPtThreshold)
            {
                //hist->Fill(tauPtFrac);
                label.type = llpdnnx::LLPLabel::Type::isPrompt_TAU;
            }

            else if (nPromptMuons > 0)
            {
                label.type = llpdnnx::LLPLabel::Type::isPrompt_MU;
            }


            else if (nPromptElectrons > 0)
            {
                label.type = llpdnnx::LLPLabel::Type::isPrompt_E;
            }

            else if (hadronFlavor==5)
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
                            label.type = llpdnnx::LLPLabel::Type::isBB;
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
                        label.type = llpdnnx::LLPLabel::Type::isCC;
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
                /*
                std::cout<<"  jet: "<<ijet<<"/"<<(jetCollection->size())<<", pt="<<jet.pt()<<", eta="<<jet.eta()<<", phi="<<jet.phi()<<", d="<<(label.displacement)<<", partonFlavor="<<(partonFlavor)<<", hadronFlavour="<<(hadronFlavor)<<", nConstituents="<<(jet.numberOfDaughters())<<std::endl;
                std::cout<<"matched gen jet: pt="<<jet.genJet()->pt()<<", eta="<<jet.genJet()->eta()<<", phi="<<jet.genJet()->phi()<<std::endl;
                for (unsigned int iConst = 0; iConst < jet.numberOfDaughters(); iConst++)
                {
                    const reco::Candidate* constituent = jet.daughter(iConst);
                    unsigned int absId = std::abs(constituent->pdgId());
                    std::cout<<"  constituent: "<<iConst<<"/"<<(jet.numberOfDaughters())<<", pt="<<constituent->pt()<<", eta="<<constituent->eta()<<", phi="<<constituent->phi()<<", pdgId="<<absId<<std::endl;
                }
                */
                label.type = llpdnnx::LLPLabel::Type::isUndefined;

            }
            
            
            
            if (displacedGenVertexCollection.product())
            {
                //find matching gen jet and corresponding vertex
                float dRmin = 0.4;
                const llpdnnx::DisplacedGenVertex* matchedVertex = nullptr;
                const reco::GenJet* matchedGenJet = nullptr;
                for (const auto& vertex: *displacedGenVertexCollection)
                {
                    for(unsigned int igenJet = 0; igenJet<vertex.genJets.size();++igenJet)
                    {
                        const reco::GenJet& genJet = vertex.genJets.at(igenJet);
                        float dRGenJets = reco::deltaR(genJet,jet);
                        
                        if(dRGenJets<dRmin)
                        {
                            dRmin = dRGenJets;
                            matchedVertex = &vertex;
                            matchedGenJet = &genJet;
                        }
                    }
                }
                
                //excludes cases where the displaced vertex is located at 0,0,0
                if (matchedVertex and matchedGenJet and matchedVertex->vertex.mag2()>(DisplacedGenVertex::MIN_DISPLACEMENT*DisplacedGenVertex::MIN_DISPLACEMENT))
                {
                    label.matchedGenJetDeltaR = dRmin;
                    label.matchedGenJetPt = matchedGenJet->pt();
                    
                    int maxFlavor = 0;
                    
                    //take displacement from matched vertex, i.e. end of potential decay chain = largest displacement
                    label.displacement = std::log10(std::max<float>(matchedVertex->d3d(),DisplacedGenVertex::MIN_DISPLACEMENT));
                    label.displacement_xy = std::log10(std::max<float>(matchedVertex->dxy(),DisplacedGenVertex::MIN_DISPLACEMENT));
                    label.displacement_z = std::log10(std::max<float>(matchedVertex->dz(),DisplacedGenVertex::MIN_DISPLACEMENT));

                    //iteratively check the decay chain to detect e.g. LLP->B->C
                    
                    const llpdnnx::DisplacedGenVertex* llpVertex = nullptr;
                    const llpdnnx::DisplacedGenVertex* parentVertex = matchedVertex;
                    while (parentVertex)
                    {
                        if (not parentVertex->motherLongLivedParticle.isNull())
                        {
                            const auto &mother = *(parentVertex->motherLongLivedParticle);
                            if (maxFlavor<getHadronFlavor(mother))
                            {
                                maxFlavor = getHadronFlavor(mother);
                                llpVertex = parentVertex;
                            }
                        }
                        parentVertex = parentVertex->motherVertex.get(); //will be null if no mother
                    }
                    if (maxFlavor>10000)
                    {
                        const std::vector<llpdnnx::LLPGhostFlavour>& llpGhostFlavours = (*llpGhostFlavourInfoMap)[jet_ref].llpFlavours;
                        
                        const auto &mother = *(llpVertex->motherLongLivedParticle);
                        label.llpId = mother.pdgId();
                        label.decay_angle = angle(matchedGenJet->p4(),mother.p4());
                        label.betagamma = mother.p()/std::max<float>(mother.mass(),DisplacedGenVertex::MIN_LLP_MASS);
                        label.llp_mass = mother.mass();
                        label.llp_pt = mother.pt();

                        int nQ = 0;
                        int nMU = 0;
                        int nE = 0;
                        int nB = 0;
                        int nTau = 0;
                        
                        for (auto const& llpGhostFlavour: llpGhostFlavours)
                        {
                            auto const& ghost = llpGhostFlavour.ghost;
                            
                            //this will associate the decay vertex with a llp decay
                            
                            if (std::sqrt((ghost->vertex()-llpVertex->vertex).mag2())>DisplacedGenVertex::MIN_DISPLACEMENT)
                            {
                                continue;
                            }
                            
                            
                            double pt = ghost->pt(); //ignore soft ghosts that wont be reconstructed
                            
                            int absId = abs(ghost->pdgId());
                            if (absId==13 and pt>muonPtThreshold_) nMU+=1;
                            if (absId==11 and pt>electronPtThreshold_) nE+=1;
                            if (absId==5 and pt>bPtThreshold_) nB+=1;
                            if ((absId<5 or absId==21) and pt>quarkPtThreshold_) nQ+=1;
                            if ((absId==15) and pt>tauPtThreshold_) nTau++;
                        }
                                
                        label.type = llpdnnx::LLPLabel::Type::isLLP_RAD; //default
                        if (nTau>=1)
                        {
                            if ((nQ+nB)==0) label.type = llpdnnx::LLPLabel::Type::isLLP_TAU;
                            if ((nQ+nB)==1) label.type = llpdnnx::LLPLabel::Type::isLLP_QTAU;
                            if ((nQ+nB)>1) label.type = llpdnnx::LLPLabel::Type::isLLP_QQTAU;
                        }
                        else if (nMU==0 and nE==0 and nTau==0)
                        {
                            if (nQ==1 and nB==0) label.type = llpdnnx::LLPLabel::Type::isLLP_Q;
                            if (nQ==0 and nB==1) label.type = llpdnnx::LLPLabel::Type::isLLP_B;
                            
                            if (nQ>1 and nB==0) label.type = llpdnnx::LLPLabel::Type::isLLP_QQ;
                            if (nB>1 or (nQ==1 and nB==1)) label.type = llpdnnx::LLPLabel::Type::isLLP_BB;
                        }
                        else if (nE>=1 and nMU==0 and nTau==0)
                        {
                            if (nQ==0 and nB==0) label.type = llpdnnx::LLPLabel::Type::isLLP_E;
                            
                            if (nQ==1 and nB==0) label.type = llpdnnx::LLPLabel::Type::isLLP_QE;
                            if (nQ==0 and nB==1) label.type = llpdnnx::LLPLabel::Type::isLLP_BE;
                            
                            if (nQ>1 and nB==0) label.type = llpdnnx::LLPLabel::Type::isLLP_QQE;
                            if (nB>1 or (nQ==1 and nB==1)) label.type = llpdnnx::LLPLabel::Type::isLLP_BBE;
                        }
                        else if (nMU>=1 and nTau==0) //accept any additional number of electrons
                        {
                            if (nQ==0 and nB==0) label.type = llpdnnx::LLPLabel::Type::isLLP_MU;
                            
                            if (nQ==1 and nB==0) label.type = llpdnnx::LLPLabel::Type::isLLP_QMU;
                            if (nQ==0 and nB==1) label.type = llpdnnx::LLPLabel::Type::isLLP_BMU;
                            
                            if (nQ>1 and nB==0) label.type = llpdnnx::LLPLabel::Type::isLLP_QQMU;
                            if (nB>1 or (nQ==1 and nB==1)) label.type = llpdnnx::LLPLabel::Type::isLLP_BBMU;
                        } 
                    }
                }
            }
        }
        //std::cout<<"  jet: "<<ijet<<"/"<<(jetCollection->size())<<", pt="<<jet.pt()<<", eta="<<jet.eta()<<", phi="<<jet.phi()<<", type="<< llpdnnx::LLPLabel::typeToString(label.type)<<", df="<<(label.llpId)<<", d="<<(label.displacement)<<std::endl;
        if (jetsPerLabel_.find(label.type)==jetsPerLabel_.end()) jetsPerLabel_[label.type] = 0;
        jetsPerLabel_[label.type] += 1;


        outputLLPLabelInfo->emplace_back(label,jet_ref);
    }

    iEvent.put(std::move(outputLLPLabelInfo));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LLPLabelOldProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(LLPLabelOldProducer);

