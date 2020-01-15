#ifndef LLPReco_DataFormats_LLPGenDecayInfo_h
#define LLPReco_DataFormats_LLPGenDecayInfo_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"

namespace llpdnnx {

struct LLPGenDecayInfo {
    std::string name;
    edm::Ptr<reco::GenParticle> llp;
    edm::PtrVector<reco::GenParticle> decayProducts;
};



}

#endif
