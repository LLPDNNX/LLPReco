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
    *(++ptr) = tag_info_features.csv_trackSumJetEtRatio;
    *(++ptr) = tag_info_features.csv_trackSumJetDeltaR;
    *(++ptr) = tag_info_features.csv_vertexCategory;
    *(++ptr) = tag_info_features.csv_trackSip2dValAboveCharm;
    *(++ptr) = tag_info_features.csv_trackSip2dSigAboveCharm;
    *(++ptr) = tag_info_features.csv_trackSip3dValAboveCharm;
    *(++ptr) = tag_info_features.csv_trackSip3dSigAboveCharm;
    *(++ptr) = tag_info_features.csv_jetNSelectedTracks;
    *(++ptr) = tag_info_features.csv_jetNTracksEtaRel;
  }

  void jet4vec_tensor_filler(tensorflow::Tensor& tensor,
                             std::size_t jet_n,
                             const llpdnnx::XTagFeatures& features) {
    float* ptr = &tensor.matrix<float>()(jet_n, 0);

    // jet 4 vector variables
    const auto& jet_features = features.jet_features;
    *ptr = jet_features.pt;
    *(++ptr) = jet_features.eta;
  }


  void cpf_tensor_filler(tensorflow::Tensor& tensor,
                          std::size_t jet_n,
                          std::size_t cpf_n,
                          const llpdnnx::ChargedCandidateFeatures& cpf_features) {
    float* ptr = &tensor.tensor<float, 3>()(jet_n, cpf_n, 0);

    *ptr = cpf_features.cpf_trackEtaRel;
    *(++ptr) = cpf_features.cpf_trackPtRel;
    *(++ptr) = cpf_features.cpf_trackPPar;
    *(++ptr) = cpf_features.cpf_trackDeltaR;
    *(++ptr) = cpf_features.cpf_trackPParRatio;
    *(++ptr) = cpf_features.cpf_trackSip2dVal;
    *(++ptr) = cpf_features.cpf_trackSip2dSig;
    *(++ptr) = cpf_features.cpf_trackSip3dVal;
    *(++ptr) = cpf_features.cpf_trackSip3dSig;
    *(++ptr) = cpf_features.cpf_trackJetDistVal;
    *(++ptr) = cpf_features.cpf_ptrel;
    *(++ptr) = cpf_features.cpf_drminsv;
    *(++ptr) = cpf_features.cpf_vertex_association;
    *(++ptr) = cpf_features.cpf_fromPV;
    *(++ptr) = cpf_features.cpf_puppi_weight;
    *(++ptr) = cpf_features.cpf_track_chi2;
    *(++ptr) = cpf_features.cpf_track_ndof;
    *(++ptr) = cpf_features.cpf_track_quality;
  }

  void npf_tensor_filler(tensorflow::Tensor& tensor,
                          std::size_t jet_n,
                          std::size_t npf_n,
                          const llpdnnx::NeutralCandidateFeatures& npf_features) {
    float* ptr = &tensor.tensor<float, 3>()(jet_n, npf_n, 0);

    *ptr = npf_features.npf_ptrel;
    *(++ptr) = npf_features.npf_deltaR;
    *(++ptr) = npf_features.npf_isGamma;
    *(++ptr) = npf_features.npf_hcal_fraction;
    *(++ptr) = npf_features.npf_drminsv;
    *(++ptr) = npf_features.npf_puppi_weight;
  }

  void sv_tensor_filler(tensorflow::Tensor& tensor,
                        std::size_t jet_n,
                        std::size_t sv_n,
                        const llpdnnx::SecondaryVertexFeatures& sv_features) {
    float* ptr = &tensor.tensor<float, 3>()(jet_n, sv_n, 0);

    *ptr = sv_features.sv_pt;
    *(++ptr) = sv_features.sv_deltaR;
    *(++ptr) = sv_features.sv_mass;
    *(++ptr) = sv_features.sv_ntracks;
    *(++ptr) = sv_features.sv_chi2;
    *(++ptr) = sv_features.sv_ndof;
    *(++ptr) = sv_features.sv_dxy;
    *(++ptr) = sv_features.sv_dxysig;
    *(++ptr) = sv_features.sv_d3d;
    *(++ptr) = sv_features.sv_d3dsig;
    *(++ptr) = sv_features.sv_costhetasvpv;
    *(++ptr) = sv_features.sv_enratio;
  }
}
