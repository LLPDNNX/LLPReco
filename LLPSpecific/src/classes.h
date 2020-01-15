#include "NANOX/LLPSpecific/interface/LLPInfo.h"
#include "DataFormats/Common/interface/Wrapper.h"

namespace { 

struct dictionary_nanox_llp
{
    nanox::LLPInfo dummy1;
    edm::Wrapper<nanox::LLPInfo> dummy2;
};

}
