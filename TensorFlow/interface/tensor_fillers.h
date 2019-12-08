#ifndef LLPReco_TensorFlow_tensor_fillers_h
#define LLPReco_TensorFlow_tensor_fillers_h

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "LLPReco/DataFormats/interface/XTagInfo.h"

namespace llpdnnx {

  // Note on setting tensor values:
  // Instead of using the more convenient tensor.matrix (etc) methods,
  // we can exploit that in the following methods values are set along
  // the innermost (= last) axis. Those values are stored contiguously in
  // the memory, so it is most performant to get the pointer to the first
  // value and use pointer arithmetic to iterate through the next pointers.

  void jet_tensor_filler(tensorflow::Tensor& tensor,
                         std::size_t jet_n,
                         const llpdnnx::XTagFeatures& features);

  void jet4vec_tensor_filler(tensorflow::Tensor& tensor,
                             std::size_t jet_n,
                             const llpdnnx::XTagFeatures& features);

  void cpf_tensor_filler(tensorflow::Tensor& tensor,
                          std::size_t jet_n,
                          std::size_t cpf_n,
                          const llpdnnx::ChargedCandidateFeatures& cpf_features);

  void npf_tensor_filler(tensorflow::Tensor& tensor,
                          std::size_t jet_n,
                          std::size_t npf_n,
                          const llpdnnx::NeutralCandidateFeatures& npf_features);

  void sv_tensor_filler(tensorflow::Tensor& tensor,
                        std::size_t jet_n,
                        std::size_t sv_n,
                        const llpdnnx::SecondaryVertexFeatures& sv_features);

}

#endif