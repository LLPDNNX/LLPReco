#include "LLPReco/XTagProducer/interface/tensor_fillers.h"

namespace llpdnnx {

  // Note on setting tensor values:
  // Instead of using the more convenient tensor.matrix (etc) methods,
  // we can exploit that in the following methods values are set along
  // the innermost (= last) axis. Those values are stored contiguously in
  // the memory, so it is most performant to get the pointer to the first
  // value and use pointer arithmetic to iterate through the next pointers.

  //    JetFeatures jet_features;

  void jet_tensor_filler(tensorflow::Tensor& tensor,
                         std::size_t jet_n,
                         const llpdnnx::XTagFeatures& features) {
    float* ptr = &tensor.matrix<float>()(jet_n, 0);

    // jet variables
    const auto& jet_features = features.jet_features;
    *ptr = jet_features.pt;
    *(++ptr) = jet_features.eta;
    // number of elements in different collections
    *(++ptr) = features.cpf_features.size();
    *(++ptr) = features.npf_features.size();
    *(++ptr) = features.sv_features.size();
    // variables from ShallowTagInfo
    const auto& tag_info_features = features.tag_info_features;
    *(++ptr) = tag_info_features.trackSumJetEtRatio;
    *(++ptr) = tag_info_features.trackSumJetDeltaR;
    *(++ptr) = tag_info_features.vertexCategory;
    *(++ptr) = tag_info_features.trackSip2dValAboveCharm;
    *(++ptr) = tag_info_features.trackSip2dSigAboveCharm;
    *(++ptr) = tag_info_features.trackSip3dValAboveCharm;
    *(++ptr) = tag_info_features.trackSip3dSigAboveCharm;
    *(++ptr) = tag_info_features.jetNSelectedTracks;
    *(++ptr) = tag_info_features.jetNTracksEtaRel;
  }

  void cpf_tensor_filler(tensorflow::Tensor& tensor,
                          std::size_t jet_n,
                          std::size_t cpf_n,
                          const llpdnnx::ChargedCandidateFeatures& cpf_features) {
    float* ptr = &tensor.tensor<float, 3>()(jet_n, cpf_n, 0);

    *ptr = cpf_features.trackEtaRel;
    *(++ptr) = cpf_features.trackPtRel;
    *(++ptr) = cpf_features.trackPPar;
    *(++ptr) = cpf_features.trackDeltaR;
    *(++ptr) = cpf_features.trackPParRatio;
    *(++ptr) = cpf_features.trackSip2dVal;
    *(++ptr) = cpf_features.trackSip2dSig;
    *(++ptr) = cpf_features.trackSip3dVal;
    *(++ptr) = cpf_features.trackSip3dSig;
    *(++ptr) = cpf_features.trackJetDistVal;
    *(++ptr) = cpf_features.ptrel;
    *(++ptr) = cpf_features.drminsv;
    *(++ptr) = cpf_features.vertex_association;
    *(++ptr) = cpf_features.fromPV;
    *(++ptr) = cpf_features.puppi_weight;
    *(++ptr) = cpf_features.track_chi2;
    *(++ptr) = cpf_features.track_ndof;
    *(++ptr) = cpf_features.track_quality;
  }

  void npf_tensor_filler(tensorflow::Tensor& tensor,
                          std::size_t jet_n,
                          std::size_t npf_n,
                          const llpdnnx::NeutralCandidateFeatures& npf_features) {
    float* ptr = &tensor.tensor<float, 3>()(jet_n, npf_n, 0);

    *ptr = npf_features.ptrel;
    *(++ptr) = npf_features.deltaR;
    *(++ptr) = npf_features.isGamma;
    *(++ptr) = npf_features.hcal_fraction;
    *(++ptr) = npf_features.drminsv;
    *(++ptr) = npf_features.puppi_weight;
  }

  void sv_tensor_filler(tensorflow::Tensor& tensor,
                        std::size_t jet_n,
                        std::size_t sv_n,
                        const llpdnnx::SecondaryVertexFeatures& sv_features) {
    float* ptr = &tensor.tensor<float, 3>()(jet_n, sv_n, 0);

    *ptr = sv_features.ptrel;
    *(++ptr) = sv_features.deltaR;
    *(++ptr) = sv_features.mass;
    *(++ptr) = sv_features.ntracks;
    *(++ptr) = sv_features.chi2;
    *(++ptr) = sv_features.ndof;
    *(++ptr) = sv_features.dxy;
    *(++ptr) = sv_features.dxysig;
    *(++ptr) = sv_features.d3d;
    *(++ptr) = sv_features.d3dsig;
    *(++ptr) = sv_features.costhetasvpv;
    *(++ptr) = sv_features.enratio;
  }
}
