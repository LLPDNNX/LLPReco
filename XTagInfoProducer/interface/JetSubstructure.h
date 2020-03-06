#ifndef LLPReco_XTagInfoProducer_JetSubstructure_h
#define LLPReco_XTagInfoProducer_JetSubstructure_h

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"


#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"

#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/MassDropTagger.hh"

#include "TLorentzVector.h"

namespace llpdnnx
{

class JetSubstructure
{
    public:
        enum class ClusterType
        {
            KT, CA, AK
        };
        
    protected:
        std::vector<fastjet::PseudoJet> consituents_;
        std::vector<TLorentzVector> lorentzVectors_;
        
        double massFromConstituents_;
        
        static fastjet::JetDefinition makeJetDefinition(ClusterType type, double r);
        
        void sortLists();
        
        std::vector<TLorentzVector> vectorsToCM() const;
        
    public:
        JetSubstructure(const reco::Jet& jet);
        JetSubstructure(const fastjet::PseudoJet& jet);
        
        inline double massFromConstituents() const
        {
            return massFromConstituents_;
        }
        
        std::vector<fastjet::PseudoJet> reclusterInclusive(ClusterType type, double r = 0.4, double minPt = 1e-10) const;
        std::vector<fastjet::PseudoJet> reclusterExclusive(ClusterType type, double r, int njets) const;
        
        //default values from https://cmssdt.cern.ch/lxr/source/RecoJets/JetProducers/python/nJettinessAdder_cfi.py except R=0.8 -> 0.4
        double nSubjettiness(
            int n,
            const fastjet::contrib::AxesDefinition& axisDef = fastjet::contrib::OnePass_KT_Axes(), 
            const fastjet::contrib::MeasureDefinition& measureDef = fastjet::contrib::NormalizedMeasure(1.0, 0.4)
        ) const;
        
        //default values from https://cmssdt.cern.ch/lxr/source/RecoJets/JetProducers/python/ak8PFJets_cfi.py
        double relMassDropMass(ClusterType type = ClusterType::AK, double r = 0.4, double muCut = 0.667, double yCut = 0.08) const;
        
        //default values from https://cmssdt.cern.ch/lxr/source/RecoJets/JetProducers/python/ak8PFJets_cfi.py
        double relSoftDropMass(ClusterType type = ClusterType::AK, double r = 0.4, double zCut = 0.1, double beta = 0.) const;
        
        
        double thrust(bool boostToCM=true) const;
        
        inline size_t nConstituents() const
        {
            return consituents_.size();
        }
        
        EventShapeVariables eventShapeVariables(bool boostToCM=true) const;
};

}

#endif
