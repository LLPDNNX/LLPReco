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
        
        static void classifyTauDecay(const reco::GenParticle* tau, std::vector<const reco::Candidate*> constituents)
        {
            
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
        
        double tauPtThreshold_;
        double quarkPtThreshold_;
        double bPtThreshold_;
        double muonPtThreshold_;
        double electronPtThreshold_;
        
        int nGenTau_;
        int nLabelledTau_;
        int nPatTau_;
        int nPatNotGen_;
        int nLabelledNotPat_;
    
        edm::EDGetTokenT<edm::View<pat::Jet>> jetToken_;
        edm::EDGetTokenT<edm::View<llpdnnx::DisplacedGenVertex>> displacedGenVertexToken_;
        edm::EDGetTokenT<edm::ValueMap<llpdnnx::LLPGhostFlavourInfo>> llpFlavourInfoToken_;
        edm::EDGetTokenT<edm::View<llpdnnx::LLPGenDecayInfo>> llpGenDecayInfoToken_;
        
        edm::EDGetTokenT<edm::View<pat::Tau>> tauToken_;

        virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

    public:
        explicit LLPLabelProducer(const edm::ParameterSet&);
        ~LLPLabelProducer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
        
        

};

LLPLabelProducer::LLPLabelProducer(const edm::ParameterSet& iConfig):
    tauPtThreshold_(iConfig.getParameter<double>("tauPtThreshold")),
    quarkPtThreshold_(iConfig.getParameter<double>("quarkPtThreshold")),
    bPtThreshold_(iConfig.getParameter<double>("bPtThreshold")),
    muonPtThreshold_(iConfig.getParameter<double>("muonPtThreshold")),
    electronPtThreshold_(iConfig.getParameter<double>("electronPtThreshold")),
    jetToken_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("srcJets"))),
    displacedGenVertexToken_(consumes<edm::View<llpdnnx::DisplacedGenVertex>>(iConfig.getParameter<edm::InputTag>("srcVertices"))),
    llpFlavourInfoToken_(consumes<edm::ValueMap<llpdnnx::LLPGhostFlavourInfo>>(iConfig.getParameter<edm::InputTag>("srcFlavourInfo"))),
    llpGenDecayInfoToken_(consumes<edm::View<llpdnnx::LLPGenDecayInfo>>( iConfig.getParameter<edm::InputTag>("srcDecayInfo") )),
    tauToken_(consumes<edm::View<pat::Tau>>( iConfig.getParameter<edm::InputTag>("srcTaus") ))
{
    produces<reco::LLPLabelInfoCollection>();
    //hist = fs->make<TH1D>("ptfrac" , "ptfrac" , 100 , 0 , 1. );
    
    nGenTau_ = 0;
    nLabelledTau_ = 0;
    nPatTau_ = 0;
    nPatNotGen_ = 0;
    nLabelledNotPat_ = 0;
}


LLPLabelProducer::~LLPLabelProducer()
{
    std::cout<<"gentau="<<nGenTau_;
    std::cout<<", labelled="<<nLabelledTau_<<" ("<<(1.*nLabelledTau_/nGenTau_);
    std::cout<<"), pat="<<nPatTau_<<" ("<<(1.*nPatTau_/nGenTau_);
    std::cout<<"), fakePat="<<nPatNotGen_<<" ("<<(1.*nPatNotGen_/nPatTau_);
    std::cout<<"), correctPat="<<(nPatTau_-nPatNotGen_)<<" ("<<(1.*(nPatTau_-nPatNotGen_)/nPatTau_);
    std::cout<<"), failedPat "<<nLabelledNotPat_<<" ("<<(1.*nLabelledNotPat_/nLabelledTau_)<<")"<<std::endl;
}




// ------------ method called to produce the data  ------------
void
LLPLabelProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    edm::Handle<edm::View<pat::Jet>> jetCollection;
    iEvent.getByToken(jetToken_, jetCollection);
          
    edm::Handle<edm::View<llpdnnx::DisplacedGenVertex>> displacedGenVertexCollection;
    iEvent.getByToken(displacedGenVertexToken_, displacedGenVertexCollection);
    
    edm::Handle<edm::ValueMap<llpdnnx::LLPGhostFlavourInfo>> llpGhostFlavourInfoMap;
    iEvent.getByToken(llpFlavourInfoToken_, llpGhostFlavourInfoMap);
    
    edm::Handle<edm::View<llpdnnx::LLPGenDecayInfo>> llpGenDecayInfoCollection;
    iEvent.getByToken(llpGenDecayInfoToken_, llpGenDecayInfoCollection);
    
    edm::Handle<edm::View<pat::Tau>> tauCollection;
    iEvent.getByToken(tauToken_, tauCollection);
    
    auto outputLLPLabelInfo = std::make_unique<reco::LLPLabelInfoCollection>();
    
    
    for (std::size_t ijet = 0; ijet < jetCollection->size(); ijet++) 
    {
        const pat::Jet& jet = jetCollection->at(ijet);
        if (jet.pt()<15.) continue;
        
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
        
            

            
            
            int nPartons = jet.jetFlavourInfo().getPartons().size();
            int nbHadrons =  jet.jetFlavourInfo().getbHadrons().size();
            int ncHadrons =  jet.jetFlavourInfo().getcHadrons().size();
            
            
            std::vector<const reco::Candidate*> promptElectrons;
            std::vector<const reco::Candidate*> promptMuons;
            std::vector<const reco::Candidate*> promptPhotons;
            
            std::vector<const reco::Candidate*> nonPromptElectrons;
            std::vector<const reco::Candidate*> nonPromptMuons;
            
            std::vector<reco::CandidatePtr> taus;
            
            for (unsigned int iConst = 0; iConst < jet.genJet()->numberOfDaughters(); iConst++)
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
                    if (absId==11) promptElectrons.push_back(packedGenConstituent);
                    if (absId==13) promptMuons.push_back(packedGenConstituent);
                    if (absId==22) promptPhotons.push_back(packedGenConstituent);
                }
                //find non prompt e/mu/photon
                else 
                {
                    if (absId==11) nonPromptElectrons.push_back(packedGenConstituent);
                    if (absId==13) nonPromptMuons.push_back(packedGenConstituent);
                }
            }
            

            const reco::Candidate* promptElectronMax = nullptr;
            const reco::Candidate* promptMuonMax = nullptr;
            const reco::Candidate* tauMax = nullptr;
            
            reco::Candidate::LorentzVector promptElectronRecoP4Max(0,0,0,0);
            reco::Candidate::LorentzVector promptMuonRecoP4Max(0,0,0,0);
            reco::Candidate::LorentzVector promptPhotonRecoP4Max(0,0,0,0);
            
            reco::Candidate::LorentzVector nonPromptElectronRecoP4Max(0,0,0,0);
            reco::Candidate::LorentzVector nonPromptMuonRecoP4Max(0,0,0,0);
            
            reco::Candidate::LorentzVector tauRecoP4Max(0,0,0,0);
            
            const auto jetConstituents = jet.getJetConstituentsQuick(); //no ref!
            
            for (const auto electron: promptElectrons)
            {
                auto match = findRecoConstituent(electron,jetConstituents);
                if (match and promptElectronRecoP4Max.pt()<match->pt()) 
                {
                    promptElectronRecoP4Max = match->p4();
                    promptElectronMax = electron;
                }
            }
            
            for (const auto muon: promptMuons)
            {
                auto match = findRecoConstituent(muon,jetConstituents);
                if (match and promptMuonRecoP4Max.pt()<match->pt())
                {
                    promptMuonRecoP4Max = match->p4();
                    promptMuonMax = muon;
                }
            }
            
            for (const auto photon: promptPhotons)
            {
                auto match = findRecoConstituent(photon,jetConstituents);
                if (match and promptPhotonRecoP4Max.pt()<match->pt())
                {
                    promptPhotonRecoP4Max = match->p4();
                }
            }
            
            for (const auto electron: nonPromptElectrons)
            {
                auto match = findRecoConstituent(electron,jetConstituents);
                if (match and nonPromptElectronRecoP4Max.pt()<match->pt()) nonPromptElectronRecoP4Max = match->p4();
            }
            
            for (const auto muon: nonPromptMuons)
            {
                auto match = findRecoConstituent(muon,jetConstituents);
                if (match and nonPromptMuonRecoP4Max.pt()<match->pt()) nonPromptMuonRecoP4Max = match->p4();
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
                if (tauRecoSum.pt()/tau->pt()>0.5)
                {
                    if (tauRecoP4Max.pt()<tauRecoSum.pt())
                    {
                        tauRecoP4Max = tauRecoSum;
                        tauMax = tau.get();
                    }
                }
            }
            
            
            //assign flavour depending on found particles
            if (nbHadrons>0)
            {
                if (nonPromptElectronRecoP4Max.pt()<1e-3 and nonPromptMuonRecoP4Max.pt()<1e-3) //note: these thresholds are arbitrary small; what counts is if the lepton is reconstructed!
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
                if (nonPromptElectronRecoP4Max.pt()<1e-3 and nonPromptMuonRecoP4Max.pt()<1e-3)
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
            else if (promptElectronRecoP4Max.pt()/jet.pt()>0.5)
            {
                label.type = llpdnnx::LLPLabel::Type::isPrompt_E;
            }
            else if (promptMuonRecoP4Max.pt()/jet.pt()>0.5)
            {
                label.type = llpdnnx::LLPLabel::Type::isPrompt_MU;
            }
            else if (promptPhotonRecoP4Max.pt()/jet.pt()>0.5)
            {
                label.type = llpdnnx::LLPLabel::Type::isPrompt_PHOTON;
            }
            else if (tauRecoP4Max.pt()/jet.pt()>0.5)
            {
                label.type = llpdnnx::LLPLabel::Type::isPrompt_TAU;
                //label.tauDecay = ?
            }
            else if (std::abs(jet.partonFlavour()!=0))
            {
                int partonFlavor = std::abs(jet.partonFlavour());
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
            
            if (taus.size()>0) nGenTau_+=1;
            
            bool matchTau = false;
            for (const auto& tau: *tauCollection)
            {
                if (reco::deltaR(tau,jet)<0.4) 
                {
                    matchTau=true;
                    nPatTau_+=1;//std::cout<<" - "<<"reco tau pt="<<tau.pt()<<std::endl;
                    if (taus.size()==0) nPatNotGen_+=1;
                }
            }
            if (label.type == llpdnnx::LLPLabel::Type::isPrompt_TAU)
            {
                nLabelledTau_+=1;
                if (not matchTau) nLabelledNotPat_+=1;
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
                        
                        
                        const reco::Candidate* displacedLepton = nullptr; //select the highest pt displaced lepton; no 50% requirement here since jet can also be q+l, qq+l
                        if (promptElectronMax)
                        {
                            //const auto electronVertex = getVertex(promptElectronMax);
                            //if (electronVertex.mag2()>DisplacedGenVertex::MIN_DISPLACEMENT and (electronVertex-llpVertex->vertex).mag2()>DisplacedGenVertex::MIN_DISPLACEMENT)
                            //{
                                if (not displacedLepton or displacedLepton->pt()<promptElectronMax->pt()) displacedLepton = promptElectronMax;
                            //}
                        }
                        
                        if (promptMuonMax)
                        {
                            //const auto muonVertex = getVertex(promptMuonMax);
                            //if (muonVertex.mag2()>DisplacedGenVertex::MIN_DISPLACEMENT and (muonVertex-llpVertex->vertex).mag2()>DisplacedGenVertex::MIN_DISPLACEMENT)
                            //{
                               if (not displacedLepton or displacedLepton->pt()<promptMuonMax->pt()) displacedLepton = promptMuonMax;
                            //}
                        }
                        
                        if (tauMax)
                        {
                            //const auto tauVertex = getVertex(tauMax);
                            //if (tauVertex.mag2()>DisplacedGenVertex::MIN_DISPLACEMENT and (tauVertex-llpVertex->vertex).mag2()>DisplacedGenVertex::MIN_DISPLACEMENT)
                            //{
                                if (not displacedLepton or displacedLepton->pt()<tauMax->pt()) displacedLepton = tauMax;
                            //}
                        }
                        
                        if (label.type==llpdnnx::LLPLabel::Type::isB or label.type==llpdnnx::LLPLabel::Type::isLeptonic_B) 
                        {
                            if (not displacedLepton)
                            {
                                label.type=llpdnnx::LLPLabel::Type::isLLP_B;
                            }
                            else
                            {
                                if (std::abs(displacedLepton->pdgId())==11) label.type=llpdnnx::LLPLabel::Type::isLLP_BE;
                                else if (std::abs(displacedLepton->pdgId())==13) label.type=llpdnnx::LLPLabel::Type::isLLP_BMU;
                                else if (std::abs(displacedLepton->pdgId())==15) label.type=llpdnnx::LLPLabel::Type::isLLP_BTAU;
                                else throw std::runtime_error("LLPLabel not well defined");
                            }
                        }
                        else if (label.type==llpdnnx::LLPLabel::Type::isBB or label.type==llpdnnx::LLPLabel::Type::isCC) 
                        {
                            if (not displacedLepton)
                            {
                                label.type=llpdnnx::LLPLabel::Type::isLLP_BB;
                            }
                            else
                            {
                                if (std::abs(displacedLepton->pdgId())==11) label.type=llpdnnx::LLPLabel::Type::isLLP_BBE;
                                else if (std::abs(displacedLepton->pdgId())==13) label.type=llpdnnx::LLPLabel::Type::isLLP_BBMU;
                                else if (std::abs(displacedLepton->pdgId())==15) label.type=llpdnnx::LLPLabel::Type::isLLP_BBTAU;
                                else throw std::runtime_error("LLPLabel not well defined");
                            }
                        }
                        else
                        {
                            if (not displacedLepton)
                            {
                                if (nPartons==0) label.type=llpdnnx::LLPLabel::Type::isLLP_RAD;
                                else if (nPartons==1) label.type=llpdnnx::LLPLabel::Type::isLLP_Q;
                                else label.type=llpdnnx::LLPLabel::Type::isLLP_QQ;
                            }
                            else
                            {
                                if (nPartons==0) 
                                {
                                    if (std::abs(displacedLepton->pdgId())==11) label.type=llpdnnx::LLPLabel::Type::isLLP_E;
                                    else if (std::abs(displacedLepton->pdgId())==13) label.type=llpdnnx::LLPLabel::Type::isLLP_MU;
                                    else if (std::abs(displacedLepton->pdgId())==15) label.type=llpdnnx::LLPLabel::Type::isLLP_TAU;
                                    else throw std::runtime_error("LLPLabel not well defined");
                                }
                                else if (nPartons==1)
                                {
                                    if (std::abs(displacedLepton->pdgId())==11) label.type=llpdnnx::LLPLabel::Type::isLLP_QE;
                                    else if (std::abs(displacedLepton->pdgId())==13) label.type=llpdnnx::LLPLabel::Type::isLLP_QMU;
                                    else if (std::abs(displacedLepton->pdgId())==15) label.type=llpdnnx::LLPLabel::Type::isLLP_QTAU;
                                    else throw std::runtime_error("LLPLabel not well defined");
                                }
                                else
                                {
                                    if (std::abs(displacedLepton->pdgId())==11) label.type=llpdnnx::LLPLabel::Type::isLLP_QQE;
                                    else if (std::abs(displacedLepton->pdgId())==13) label.type=llpdnnx::LLPLabel::Type::isLLP_QQMU;
                                    else if (std::abs(displacedLepton->pdgId())==15) label.type=llpdnnx::LLPLabel::Type::isLLP_QQTAU;
                                    else throw std::runtime_error("LLPLabel not well defined");
                                }
                            }
                        }
                    }
                }
            }
            if (label.type!=llpdnnx::LLPLabel::Type::isPU)
            std::cout<<"  jet: "<<ijet<<"/"<<(jetCollection->size())<<", pt="<<jet.pt()<<", eta="<<jet.eta()<<", type="<< llpdnnx::LLPLabel::typeToString(label.type)<<std::endl;
        }


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
