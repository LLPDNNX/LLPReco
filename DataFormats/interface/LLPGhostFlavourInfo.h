#ifndef LLPReco_DataFormats_LLPGhostFlavourInfo_h
#define LLPReco_DataFormats_LLPGhostFlavourInfo_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "LLPReco/DataFormats/interface/LLPGenDecayInfo.h"

namespace llpdnnx {

struct LLPGhostFlavour
{
    edm::Ptr<llpdnnx::LLPGenDecayInfo> decay;
    edm::Ptr<reco::GenParticle> ghost;
    
    LLPGhostFlavour(
        edm::Ptr<llpdnnx::LLPGenDecayInfo> decay,
        edm::Ptr<reco::GenParticle> ghost
    ):
        decay(decay),
        ghost(ghost)
    {
    }
};

struct LLPGhostFlavourInfo
{
    std::vector<LLPGhostFlavour> llpFlavours;
};

}

#endif
