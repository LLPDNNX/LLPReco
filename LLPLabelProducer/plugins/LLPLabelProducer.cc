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
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "LLPReco/DataFormats/interface/LLPGhostFlavourInfo.h"

#include "TH1.h"


using llpdnnx::DisplacedGenVertex;
using llpdnnx::DisplacedGenVertexCollection;

class LLPLabelProducer:
    public edm::stream::EDProducer<>
    
{
    private:
        struct CandidateHash
        {
            long operator() (const reco::CandidatePtr& cand) const 
            {
                return cand.id().id() * 100000 + cand.key();
            }
        };
        
    
        struct LLPDecay
        {
            const reco::GenParticle* llp;
            std::vector<const reco::GenParticle*> decayProducts;
            
            LLPDecay():
                llp(nullptr)
            {
            }
        };
        
        static void gatherVisibleDecayProducts(std::vector<const reco::Candidate*>& list, const reco::Candidate* genParticle)
        {
            if (genParticle->numberOfDaughters()==0) 
            {
                int absId = std::abs(genParticle->pdgId());
                if (absId!=12 and absId!=14 and absId!=16 and absId>6 and absId<100000)
                {
                    list.push_back(genParticle);
                }
            }
            else
            {
                for (size_t d=0; d< genParticle->numberOfDaughters(); ++d)
                {
                    gatherVisibleDecayProducts(list, genParticle->daughter(d));
                }
            }
        }
        
        static reco::Candidate::Point getVertex(const reco::Candidate* p)
        {
            const auto genParticle = dynamic_cast<const reco::GenParticle*>(p);
            if (genParticle) return genParticle->vertex();
            const auto packedGenParticle = dynamic_cast<const pat::PackedGenParticle*>(p);
            if (packedGenParticle) return packedGenParticle->vertex();
            return reco::Candidate::Point(0,0,0);             
        }
        
        static const reco::Candidate* findRecoConstituent(const reco::Candidate* genParticle, const std::vector<const reco::Candidate*>& constituents)
        {
            double minDR2 = 1000000;
            const reco::Candidate* match = nullptr;
            for (const auto& constituent: constituents)
            {
                double ptRel = std::fabs(constituent->pt()/genParticle->pt()-1.);
                if (ptRel>0.5) continue; //prevent random match
                double dEta = std::fabs(genParticle->eta()-constituent->eta());
                if (dEta>0.05) continue;
                double dPhi = std::fabs(reco::deltaPhi(*genParticle,*constituent));
                if (dPhi>0.1) continue;
                double dR2 = dEta*dEta+0.5*dPhi*dPhi; //allow for more spread in phi
                if (minDR2>dR2)
                {
                    minDR2 = dR2;
                    match = constituent;
                }
            }
            return match;
        }
        
        static double calcFraction(const reco::Candidate::LorentzVector& p, const reco::Candidate::LorentzVector& base)
        {
            return p.Vect().Dot(base.Vect())/base.Vect().mag2();
        }
        
        static void printDecayChain(std::vector<const reco::Candidate*>& finals, const reco::Candidate* p, int indent=1)
        {
            for (int i=0; i < indent; ++i) std::cout<<" - ";
            std::cout<<"id="<<(p->pdgId())<<", pt="<<(p->pt())<<std::endl;
            if (p->numberOfDaughters()==0) finals.push_back(p);
            for (size_t d=0; d< p->numberOfDaughters(); ++d)
            {
                printDecayChain(finals,p->daughter(d), indent+1);
            }
        }

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

        
        std::unordered_map<llpdnnx::LLPLabel::Type,int> jetsPerLabel_;
    
        edm::EDGetTokenT<edm::View<pat::Jet>> jetToken_;
        edm::EDGetTokenT<edm::View<llpdnnx::DisplacedGenVertex>> displacedGenVertexToken_;
        edm::EDGetTokenT<edm::View<llpdnnx::LLPGenDecayInfo>> llpGenDecayInfoToken_;
        

        virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

    public:
        explicit LLPLabelProducer(const edm::ParameterSet&);
        ~LLPLabelProducer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
        
        

};

LLPLabelProducer::LLPLabelProducer(const edm::ParameterSet& iConfig):
    jetToken_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("srcJets"))),
    displacedGenVertexToken_(consumes<edm::View<llpdnnx::DisplacedGenVertex>>(iConfig.getParameter<edm::InputTag>("srcVertices"))),
    llpGenDecayInfoToken_(consumes<edm::View<llpdnnx::LLPGenDecayInfo>>( iConfig.getParameter<edm::InputTag>("srcDecayInfo") ))
{
    produces<reco::LLPLabelInfoCollection>();
}


LLPLabelProducer::~LLPLabelProducer()
{
    for (const auto& labelCountPair: jetsPerLabel_)
    {
        std::cout<<llpdnnx::LLPLabel::typeToString(labelCountPair.first)<<": "<<labelCountPair.second<<std::endl;
    }
}




// ------------ method called to produce the data  ------------
void
LLPLabelProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    edm::Handle<edm::View<pat::Jet>> jetCollection;
    iEvent.getByToken(jetToken_, jetCollection);
          
    edm::Handle<edm::View<llpdnnx::DisplacedGenVertex>> displacedGenVertexCollection;
    iEvent.getByToken(displacedGenVertexToken_, displacedGenVertexCollection);
    
    edm::Handle<edm::View<llpdnnx::LLPGenDecayInfo>> llpGenDecayInfoCollection;
    iEvent.getByToken(llpGenDecayInfoToken_, llpGenDecayInfoCollection);
    
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
            const auto& genJet = jet.genJet();
            
            label.partonFlavor = std::abs(jet.partonFlavour());
            label.hadronFlavor = std::abs(jet.hadronFlavour());
            label.matchedGenJetDeltaR = reco::deltaR(jet.p4(),jet.genJet()->p4());
            label.matchedGenJetPt = jet.genJet()->pt();
        
            //std::cout<<"jet "<<ijet<<": pt="<<jet.pt()<<", genpt="<<genJet->pt()<<std::endl;
        
            int nbHadrons =  jet.jetFlavourInfo().getbHadrons().size();
            int ncHadrons =  jet.jetFlavourInfo().getcHadrons().size();
            
            
            std::vector<reco::CandidatePtr> promptElectrons;
            std::vector<reco::CandidatePtr> promptMuons;
            std::vector<reco::CandidatePtr> promptPhotons;
            
            std::vector<reco::CandidatePtr> nonPromptElectrons;
            std::vector<reco::CandidatePtr> nonPromptMuons;
            
            std::vector<reco::CandidatePtr> taus;
            
            for (unsigned int iConst = 0; iConst < genJet->numberOfDaughters(); iConst++)
            {
                const pat::PackedGenParticle* packedGenConstituent = dynamic_cast<const pat::PackedGenParticle*>(genJet->daughter(iConst));
                if (not packedGenConstituent) continue;
                if (not (packedGenConstituent->mother(0))) continue;
                
                unsigned int absId = std::abs(packedGenConstituent->pdgId());
                
                //find tau decay products
                if (packedGenConstituent->statusFlags().isHardProcessTauDecayProduct())
                {
                    
                    auto motherRef = packedGenConstituent->motherRef();
                    while (motherRef.isNonnull() and motherRef->numberOfMothers()>0 and std::fabs(motherRef->pdgId())!=15)
                    {
                        motherRef = motherRef->motherRef(0);
                    }
                    if (motherRef.isNonnull()) 
                    {
                        reco::CandidatePtr candidatePtr(motherRef.id(),motherRef.get(),motherRef.key());
                        if (std::find(taus.begin(),taus.end(),candidatePtr)==taus.end())
                        {
                            taus.push_back(candidatePtr);
                        }
                    }
                }
                //find prompt e/mu/photon
                else if (packedGenConstituent->isPromptFinalState())
                {
                    if (absId==11) promptElectrons.push_back(genJet->daughterPtr(iConst));
                    if (absId==13) promptMuons.push_back(genJet->daughterPtr(iConst));
                    if (absId==22) promptPhotons.push_back(genJet->daughterPtr(iConst));
                }
                //find non prompt e/mu/photon
                else 
                {
                    if (absId==11) nonPromptElectrons.push_back(genJet->daughterPtr(iConst));
                    if (absId==13) nonPromptMuons.push_back(genJet->daughterPtr(iConst));
                }
            }
            

            int promptRecoMaxId = -1;
            reco::Candidate::LorentzVector promptRecoMax(0,0,0,0);
            reco::Candidate::LorentzVector nonPromptLeptonRecoMax(0,0,0,0);
            
            const auto jetConstituents = jet.getJetConstituentsQuick(); //no ref!
            
            for (const auto electron: promptElectrons)
            {
                auto match = findRecoConstituent(electron.get(),jetConstituents);
                if (match) 
                {
                    if (calcFraction(promptRecoMax,jet.p4())<calcFraction(match->p4(),jet.p4()))
                    {
                        promptRecoMaxId = 11;
                        promptRecoMax = match->p4();
                    }
                }
            }
            
            for (const auto muon: promptMuons)
            {
                auto match = findRecoConstituent(muon.get(),jetConstituents);
                if (match)
                {
                    if (calcFraction(promptRecoMax,jet.p4())<calcFraction(match->p4(),jet.p4()))
                    {
                        promptRecoMaxId = 13;
                        promptRecoMax = match->p4();
                    }
                }
            }
            
            for (const auto photon: promptPhotons)
            {
                auto match = findRecoConstituent(photon.get(),jetConstituents);
                if (match)
                {
                    if (calcFraction(promptRecoMax,jet.p4())<calcFraction(match->p4(),jet.p4()))
                    {
                        promptRecoMaxId = 22;
                        promptRecoMax = match->p4();
                    }
                }
            }
            
            for (const auto electron: nonPromptElectrons)
            {
                auto match = findRecoConstituent(electron.get(),jetConstituents);
                if (match and calcFraction(nonPromptLeptonRecoMax,jet.p4())<calcFraction(nonPromptLeptonRecoMax,match->p4())) nonPromptLeptonRecoMax = match->p4();
            }
            
            for (const auto muon: nonPromptMuons)
            {
                auto match = findRecoConstituent(muon.get(),jetConstituents);
                if (match and calcFraction(nonPromptLeptonRecoMax,jet.p4())<calcFraction(nonPromptLeptonRecoMax,match->p4())) nonPromptLeptonRecoMax = match->p4();
            }
            
            for (const auto tau: taus)
            {
                reco::Candidate::LorentzVector tauRecoSum(0,0,0,0);
                std::vector<const reco::Candidate*> decayProducts;
                gatherVisibleDecayProducts(decayProducts, tau.get());
                for (const auto decayProduct: decayProducts)
                {
                    auto match = findRecoConstituent(decayProduct,jetConstituents);
                    if (match) tauRecoSum+=match->p4(); 
                }
                if (calcFraction(tauRecoSum,tau->p4())>0.5) //at least 50% of momentum reconstructed
                {
                    if (calcFraction(promptRecoMax,jet.p4())<calcFraction(tauRecoSum,jet.p4()))
                    {
                        promptRecoMaxId = 15;
                        promptRecoMax = tauRecoSum;
                    }
                    
                    //classify tau decay
                    int nH = 0;
                    int nPi0 = 0;
                    int nE = 0;
                    int nMu = 0;
                    for (size_t idaughter = 0; idaughter < tau->numberOfDaughters(); ++idaughter)
                    {
                        int daughterAbsId = std::abs( tau->daughter(idaughter)->pdgId());
                        if (daughterAbsId==11) nE+=1;
                        else if (daughterAbsId==13) nMu+=1;
                        else if (daughterAbsId==111) nPi0+=1;
                        else if (daughterAbsId>111 and daughterAbsId<500) nH+=1;
                    }
                    if (nE>0) label.tauDecay = llpdnnx::LLPLabel::TauDecay::E;
                    else if (nMu>0) label.tauDecay = llpdnnx::LLPLabel::TauDecay::MU;
                    else if (nH==1) 
                    {
                        if (nPi0==0) label.tauDecay = llpdnnx::LLPLabel::TauDecay::H;
                        else if (nPi0==1) label.tauDecay = llpdnnx::LLPLabel::TauDecay::H_1PI0;
                        else label.tauDecay = llpdnnx::LLPLabel::TauDecay::H_XPI0;
                    }   
                    else if (nH>1) 
                    {
                        if (nPi0==0) label.tauDecay = llpdnnx::LLPLabel::TauDecay::HHH;
                        else label.tauDecay = llpdnnx::LLPLabel::TauDecay::HHH_XPI0;
                    }   
                    else
                    {
                        std::cout<<"WARNING: Tau decay cannot be characterised"<<std::endl;
                    }
                }
                else
                {
                    label.tauDecay = llpdnnx::LLPLabel::TauDecay::INVISIBLE;
                }
            }
            
            
            //assign flavour depending on found particles
            if (nbHadrons>0)
            {
                if (nonPromptLeptonRecoMax.pt()<1e-3) //note: these thresholds are arbitrary small; what counts is if the lepton is reconstructed!
                {
                    if (nbHadrons>1)
                    {
                        label.type = llpdnnx::LLPLabel::Type::isBB;
                    }
                    else
                    {
                        label.type = llpdnnx::LLPLabel::Type::isB;
                    }
                }
                else
                {
                    label.type = llpdnnx::LLPLabel::Type::isLeptonic_B;
                }
            }
            else if (nbHadrons==0 and ncHadrons>0)
            {
                if (nonPromptLeptonRecoMax.pt()<1e-3)
                {
                    if (ncHadrons>1)
                    {
                        label.type = llpdnnx::LLPLabel::Type::isCC;
                    }
                    else
                    {
                        label.type = llpdnnx::LLPLabel::Type::isC;
                    }
                }
                else
                {
                    label.type = llpdnnx::LLPLabel::Type::isLeptonic_C;
                }
            }
            else if (calcFraction(promptRecoMax,jet.p4())>0.5)
            {
                if (promptRecoMaxId==11) label.type = llpdnnx::LLPLabel::Type::isPrompt_E;
                else if (promptRecoMaxId==13) label.type = llpdnnx::LLPLabel::Type::isPrompt_MU;
                else if (promptRecoMaxId==15) label.type = llpdnnx::LLPLabel::Type::isPrompt_TAU;
                else if (promptRecoMaxId==22) label.type = llpdnnx::LLPLabel::Type::isPrompt_PHOTON;
            }
            else if (std::abs(jet.partonFlavour()!=0))
            {
                int partonFlavor = std::abs(jet.partonFlavour());
                if (partonFlavor==5) 
                {
                    if (nonPromptLeptonRecoMax.pt()<1e-3) label.type = llpdnnx::LLPLabel::Type::isLeptonic_B;
                    else label.type = llpdnnx::LLPLabel::Type::isB;
                }
                if (partonFlavor==4)
                {
                    if (nonPromptLeptonRecoMax.pt()<1e-3) label.type = llpdnnx::LLPLabel::Type::isLeptonic_C;
                    else label.type = llpdnnx::LLPLabel::Type::isC;
                }
                if (partonFlavor==3) label.type = llpdnnx::LLPLabel::Type::isS;
                if (partonFlavor==2 or partonFlavor==1) label.type = llpdnnx::LLPLabel::Type::isUD;
                if (partonFlavor==21) label.type = llpdnnx::LLPLabel::Type::isG;
            }
            else //no hadrons and parton flavour=0 => jet cannot be classified!
            {
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
                        
                        const auto &mother = *(llpVertex->motherLongLivedParticle);
                        label.llpId = mother.pdgId();
                        label.decay_angle = angle(matchedGenJet->p4(),mother.p4());
                        label.betagamma = mother.p()/std::max<float>(mother.mass(),DisplacedGenVertex::MIN_LLP_MASS);
                        label.llp_mass = mother.mass();
                        label.llp_pt = mother.pt();
                        
                        
                        int nQuarks = 0;
                        for (const auto& llpGenDecayInfo: *llpGenDecayInfoCollection)
                        {
                            if ((llpGenDecayInfo.decayProducts[0]->vertex()-llpVertex->vertex).mag2()>DisplacedGenVertex::MIN_DISPLACEMENT) continue;
                            for (const auto& particle: llpGenDecayInfo.decayProducts)
                            {
                                if (std::abs(particle->pdgId())<5 and reco::deltaR(*particle,jet)<0.4) nQuarks+=1;
                            }
                        }
                            
                        if (label.type==llpdnnx::LLPLabel::Type::isB or label.type==llpdnnx::LLPLabel::Type::isLeptonic_B) 
                        {
                            if (calcFraction(promptRecoMax,jet.p4())<0.1) //arbitrary threshold; reconstructed candidates count!
                            {
                                label.type=llpdnnx::LLPLabel::Type::isLLP_B;
                            }
                            else
                            {
                                if (promptRecoMaxId==11) label.type=llpdnnx::LLPLabel::Type::isLLP_BE;
                                else if (promptRecoMaxId==13) label.type=llpdnnx::LLPLabel::Type::isLLP_BMU;
                                else if (promptRecoMaxId==15) label.type=llpdnnx::LLPLabel::Type::isLLP_BTAU;
                                else if (promptRecoMaxId==22) label.type=llpdnnx::LLPLabel::Type::isLLP_BPHOTON;
                                else throw std::runtime_error("LLPLabel not well defined");
                            }
                        }
                        
                        else if (label.type==llpdnnx::LLPLabel::Type::isBB or label.type==llpdnnx::LLPLabel::Type::isCC) 
                        {
                            if (calcFraction(promptRecoMax,jet.p4())<0.1) //arbitrary threshold; reconstructed candidates count!
                            {
                                label.type=llpdnnx::LLPLabel::Type::isLLP_BB;
                            }
                            else
                            {
                                if (promptRecoMaxId==11) label.type=llpdnnx::LLPLabel::Type::isLLP_BBE;
                                else if (promptRecoMaxId==13) label.type=llpdnnx::LLPLabel::Type::isLLP_BBMU;
                                else if (promptRecoMaxId==15) label.type=llpdnnx::LLPLabel::Type::isLLP_BBTAU;
                                else if (promptRecoMaxId==22) label.type=llpdnnx::LLPLabel::Type::isLLP_BBPHOTON;
                                else throw std::runtime_error("LLPLabel not well defined");
                            }
                        }
                        else
                        {
                            if (calcFraction(promptRecoMax,jet.p4())<0.1) //arbitrary threshold; reconstructed candidates count!
                            {
                                if (nQuarks==0) label.type=llpdnnx::LLPLabel::Type::isLLP_RAD;
                                else if (nQuarks==1) label.type=llpdnnx::LLPLabel::Type::isLLP_Q;
                                else label.type=llpdnnx::LLPLabel::Type::isLLP_QQ;
                            }
                            else
                            {
                                if (nQuarks==0) 
                                {
                                    if (promptRecoMaxId==11) label.type=llpdnnx::LLPLabel::Type::isLLP_E;
                                    else if (promptRecoMaxId==13) label.type=llpdnnx::LLPLabel::Type::isLLP_MU;
                                    else if (promptRecoMaxId==15) label.type=llpdnnx::LLPLabel::Type::isLLP_TAU;
                                    else if (promptRecoMaxId==22) label.type=llpdnnx::LLPLabel::Type::isLLP_PHOTON;
                                    else throw std::runtime_error("LLPLabel not well defined");
                                }
                                else if (nQuarks==1)
                                {
                                    if (promptRecoMaxId==11) label.type=llpdnnx::LLPLabel::Type::isLLP_QE;
                                    else if (promptRecoMaxId==13) label.type=llpdnnx::LLPLabel::Type::isLLP_QMU;
                                    else if (promptRecoMaxId==15) label.type=llpdnnx::LLPLabel::Type::isLLP_QTAU;
                                    else if (promptRecoMaxId==22) label.type=llpdnnx::LLPLabel::Type::isLLP_QPHOTON;
                                    else throw std::runtime_error("LLPLabel not well defined");
                                }
                                else
                                {
                                    if (promptRecoMaxId==11) label.type=llpdnnx::LLPLabel::Type::isLLP_QQE;
                                    else if (promptRecoMaxId==13) label.type=llpdnnx::LLPLabel::Type::isLLP_QQMU;
                                    else if (promptRecoMaxId==15) label.type=llpdnnx::LLPLabel::Type::isLLP_QQTAU;
                                    else if (promptRecoMaxId==22) label.type=llpdnnx::LLPLabel::Type::isLLP_QQPHOTON;
                                    else throw std::runtime_error("LLPLabel not well defined");
                                }
                            }
                        }
                    }
                }
            }
        }
        if (jetsPerLabel_.find(label.type)==jetsPerLabel_.end()) jetsPerLabel_[label.type] = 0;
        jetsPerLabel_[label.type] += 1;

        outputLLPLabelInfo->emplace_back(label,jet_ref);
    }

    iEvent.put(std::move(outputLLPLabelInfo));
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
