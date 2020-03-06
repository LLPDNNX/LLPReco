#include "LLPReco/XTagInfoProducer/interface/JetSubstructure.h"

namespace llpdnnx
{

JetSubstructure::JetSubstructure(const reco::Jet& jet)
{
    if (jet.numberOfDaughters()==0) throw cms::Exception("Jets has no constituents!");
    TLorentzVector jetVectorFromConsituents(0,0,0,0);
    for(unsigned int iconstituent = 0; iconstituent < jet.numberOfDaughters(); ++iconstituent)
    {
        const reco::Candidate* constituent = jet.daughter(iconstituent);
        if((constituent->energy() < 1e-10) or (constituent->mass()<1e-10))
        {
            continue;
        }
        consituents_.emplace_back(constituent->px(),constituent->py(),constituent->pz(),constituent->energy());
        lorentzVectors_.emplace_back(constituent->px(),constituent->py(),constituent->pz(),constituent->energy());
        jetVectorFromConsituents+=lorentzVectors_.back();
    }
    
    massFromConstituents_ = std::max(1e-10,jetVectorFromConsituents.M());
    
    sortLists();
}

JetSubstructure::JetSubstructure(const fastjet::PseudoJet& jet)
{
    if (jet.constituents().size()==0) throw cms::Exception("Jets has no constituents!");
    TLorentzVector jetVectorFromConsituents(0,0,0,0);
    for(auto const& constituent: jet.constituents())
    {
        if((constituent.e()<1e-10) or (constituent.m()<1e-10))
        {
            continue;
        }
        consituents_.emplace_back(constituent);
        lorentzVectors_.emplace_back(constituent.px(),constituent.py(),constituent.pz(),constituent.e());
        jetVectorFromConsituents+=lorentzVectors_.back();
    }
    
    massFromConstituents_ = std::max(1e-10,jetVectorFromConsituents.M());
    
    sortLists();
}

void JetSubstructure::sortLists()
{
    std::sort(consituents_.begin(),consituents_.end(),[](const fastjet::PseudoJet& a, const fastjet::PseudoJet& b) 
    { 
        return a.modp2() > b.modp2();
    });
    std::sort(lorentzVectors_.begin(),lorentzVectors_.end(),[](const TLorentzVector& a, const TLorentzVector& b) 
    { 
        return a.Mag2() > b.Mag2();
    });
}

std::vector<TLorentzVector> JetSubstructure::vectorsToCM() const
{
    TLorentzVector cm(0,0,0,0);
    for (auto const& lorentzVector: lorentzVectors_)
    {
        cm += lorentzVector;
    }
    
    auto boostVector = cm.BoostVector();
    
    std::vector<TLorentzVector> boostedVectors;
    for (auto const& lorentzVector: lorentzVectors_)
    { 
        TLorentzVector boostedVector = lorentzVector;
        boostedVector.Boost(-boostVector);
        boostedVectors.push_back(boostedVector);
    }
    return boostedVectors;
}

fastjet::JetDefinition JetSubstructure::makeJetDefinition(ClusterType type, double r)
{
    switch (type)
    {
        case JetSubstructure::ClusterType::KT:
            return fastjet::JetDefinition(fastjet::kt_algorithm, r);
        case JetSubstructure::ClusterType::CA:
            return fastjet::JetDefinition(fastjet::cambridge_algorithm, r);
        case JetSubstructure::ClusterType::AK:
            return fastjet::JetDefinition(fastjet::antikt_algorithm, r);
    }
    return fastjet::JetDefinition(fastjet::antikt_algorithm, r);
}

std::vector<fastjet::PseudoJet> JetSubstructure::reclusterInclusive(
    JetSubstructure::ClusterType type, 
    double r, 
    double minPt
) const
{
    fastjet::ClusterSequence clusterSequence(consituents_, makeJetDefinition(type,r));
    std::vector<fastjet::PseudoJet> reclusteredJets = fastjet::sorted_by_pt(
        clusterSequence.inclusive_jets(minPt)
    );
    
    return reclusteredJets;
}

std::vector<fastjet::PseudoJet> JetSubstructure::reclusterExclusive(
    JetSubstructure::ClusterType type, 
    double r, 
    int njets
) const
{
    fastjet::ClusterSequence clusterSequence(consituents_, makeJetDefinition(type,r));
    std::vector<fastjet::PseudoJet> reclusteredJets = fastjet::sorted_by_pt(
        clusterSequence.exclusive_jets(njets)
    );
    
    return reclusteredJets;
}

double JetSubstructure::nSubjettiness(
    int n,
    const fastjet::contrib::AxesDefinition& axisDef, 
    const fastjet::contrib::MeasureDefinition& measureDef
) const
{
    fastjet::contrib::Njettiness njettiness(axisDef,measureDef);
    return njettiness.getTau(n, consituents_);
}

double JetSubstructure::thrust(bool boostToCM) const
{
    std::vector<TLorentzVector> particleVectors;

    if (lorentzVectors_.size() < 3) return 0;
    
    if (boostToCM)
    {
        particleVectors = vectorsToCM();
    }
    else
    {
        particleVectors = lorentzVectors_;
    }
    
    unsigned int n = std::min<unsigned int>(particleVectors.size(),4); //3 or 4
    
    std::vector<TVector3> tvec;
    std::vector<double> tval;
    
    //all permutations of n leading jets
    for (int i = 0 ; i < std::pow(2, n-1); ++i) 
    {
        // Create an initial guess from the n leading jets
        TVector3 thrustAxis(0,0,0);
        int sign = i;
        for (unsigned int k = 0 ; k < n ; ++k)
        {
            if ((sign % 2) == 1)
            {
                thrustAxis += particleVectors[k].Vect();
            }
            else
            {
                thrustAxis -= particleVectors[k].Vect();
            }
            sign /= 2;
        }
        thrustAxis=thrustAxis.Unit();

        // Iterate
        double diff=999.;
        while (diff>1e-5) 
        {
            TVector3 diffVec(0,0,0);
            for (unsigned int k = 0 ; k<particleVectors.size() ; k++)
            {
                if (thrustAxis.Dot(particleVectors[k].Vect())>0)
                {
                    diffVec+=particleVectors[k].Vect();
                }
                else
                {
                    diffVec-=particleVectors[k].Vect();
                }
            }
            diff=(thrustAxis-diffVec.Unit()).Mag();
            thrustAxis=diffVec.Unit();
        }

        // Calculate the thrust value for the vector we found
        double thrust = 0.;
        double sum = 0.;
        for (unsigned int k=0 ; k<particleVectors.size() ; k++)
        {
            thrust+=fabs(thrustAxis.Dot(particleVectors[k].Vect()));
            sum+=particleVectors[k].Vect().Mag();
        }
        thrust /= sum;

        tval.push_back(thrust);
        tvec.push_back(thrustAxis);
    }

    // Pick the solution with the largest thrust
    double maxThrust = 0.;
    for (unsigned int i=0 ; i<tvec.size() ; i++)
    {
        if (tval[i] > maxThrust)
        {
            maxThrust = tval[i];
        }
    }
    return maxThrust;
}

double JetSubstructure::relMassDropMass(ClusterType type, double r, double muCut, double yCut) const
{
    if (massFromConstituents_<1e-10) return 0;
    if (consituents_.size()<2) return 0;
    
    fastjet::ClusterSequence clusterSequence(consituents_, makeJetDefinition(type,r));
    std::vector<fastjet::PseudoJet> reclusteredJets = fastjet::sorted_by_pt(
        clusterSequence.inclusive_jets(1e-10)
    );
    if (reclusteredJets.size()==0) throw cms::Exception("No jets found after reclustering!");
    if (reclusteredJets[0].constituents().size()<2) return 0;
    fastjet::MassDropTagger massDropTagger(muCut, yCut);
    fastjet::PseudoJet taggedJet = massDropTagger.result(reclusteredJets[0]);
    
    double massDropMass = taggedJet.m();
    if (massDropMass<1e-10) return 0;
    
    return massDropMass/massFromConstituents_;
}

double JetSubstructure::relSoftDropMass(ClusterType type, double r, double zCut, double beta) const
{
    if (massFromConstituents_<1e-10) return 0;
    if (consituents_.size()<2) return 0;
    
    fastjet::ClusterSequence clusterSequence(consituents_, makeJetDefinition(type,r));
    std::vector<fastjet::PseudoJet> reclusteredJets = fastjet::sorted_by_pt(
        clusterSequence.inclusive_jets(1e-10)
    );
    if (reclusteredJets.size()==0) throw cms::Exception("No jets found after reclustering!");
    if (reclusteredJets[0].constituents ().size()<2) return 0;
    
    fastjet::contrib::SoftDrop softDrop(beta, zCut, r );
    fastjet::PseudoJet softDropJet = softDrop.result(reclusteredJets[0]);
    
    double softDropMass = softDropJet.m();
    if (softDropMass<1e-10) return 0;
    
    return softDropMass/massFromConstituents_;
}

EventShapeVariables JetSubstructure::eventShapeVariables(bool boostToCM) const
{
    std::vector<math::XYZVector> vectors;
    
    if (boostToCM)
    {
        for (auto const& lorentzVector: vectorsToCM())
        {
            vectors.emplace_back(lorentzVector.Px(),lorentzVector.Py(),lorentzVector.Pz());
        }
    }
    else
    {
        for (auto const& lorentzVector: lorentzVectors_)
        {
            vectors.emplace_back(lorentzVector.Px(),lorentzVector.Py(),lorentzVector.Pz());
        }
    }
    return EventShapeVariables(vectors);
}

}

