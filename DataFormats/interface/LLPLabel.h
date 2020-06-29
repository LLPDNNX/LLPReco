#ifndef LLPReco_DataFormats_LLPLabel_h
#define LLPReco_DataFormats_LLPLabel_h

#include "DataFormats/BTauReco/interface/FeaturesTagInfo.h"
#include "LLPReco/DataFormats/interface/DisplacedGenVertex.h"

namespace llpdnnx {

class LLPLabel {

    public:

    
        enum class Type {
            isPU,
            isB,
            isBB,
            isGBB,
            isLeptonic_B,
            isLeptonic_C,
            isB_MU,
            isB_E,
            isC_MU,
            isC_E,
            isC,
            isCC,
            isGCC,
            isS,
            isS_MU,
            isS_E,
            isUD,
            isUD_MU,
            isUD_E,
            isG,
            isG_MU,
            isG_E,
            
            isLLP_RAD, //no flavour match (likely from wide angle radiation)
            isLLP_MU, //prompt lepton
            isLLP_E, //prompt lepton
            isLLP_Q, //single light quark
            isLLP_QMU, //single light quark + prompt lepton
            isLLP_QE, //single light quark + prompt lepton
            isLLP_QQ, //double light quark
            isLLP_QQMU, //double light quark + prompt lepton
            isLLP_QQE, //double light quark + prompt lepton
            isLLP_B, //single b/c quark
            isLLP_BMU, //single b/c quark + prompt lepton
            isLLP_BE, //single b/c quark + prompt lepton
            isLLP_BB, //double b/c quark
            isLLP_BBMU, //double b/c quark + prompt lepton
            isLLP_BBE, //double b/c quark + prompt lepton
            isLLP_TAU,
            isLLP_QTAU,
            isLLP_QQTAU,
            isUndefined
            
        };
        
        Type type;
        
        int jetIdx;
        int partonFlavor;
        int hadronFlavor;
        int llpId;
        float llp_mass;
        float llp_pt;

        float displacement;
        float displacement_xy;
        float displacement_z;
        float decay_angle;
        float betagamma;
        
        float matchedGenJetDeltaR;
        float matchedGenJetPt;
        
        LLPLabel():
            type(Type::isUndefined),
            partonFlavor(0),
            hadronFlavor(0),
            llpId(0),
            llp_mass(DisplacedGenVertex::MIN_LLP_MASS),
            llp_pt(0),
            displacement(std::log10(DisplacedGenVertex::MIN_DISPLACEMENT)), //log10(x/1cm)
            displacement_xy(std::log10(DisplacedGenVertex::MIN_DISPLACEMENT)), //log10(x/1cm)
            displacement_z(std::log10(DisplacedGenVertex::MIN_DISPLACEMENT)), //log10(x/1cm)
            decay_angle(0),
            betagamma(0),
            matchedGenJetDeltaR(-1),
            matchedGenJetPt(-1)
        {
        }
        
        inline static const std::string typeToString(const Type& type)
        {
            switch (type)
            {
                case Type::isPU:
                    return "isPU";
                case Type::isB:
                    return "isB";
                case Type::isBB:
                    return "isBB";
                case Type::isGBB:
                    return "isGBB";
                case Type::isLeptonic_B:
                    return "isLeptonic_B";
                case Type::isLeptonic_C:
                    return "isLeptonic_C";
                case Type::isB_MU:
                    return "isB_MU";
                case Type::isB_E:
                    return "isB_E";
                case Type::isC_MU:
                    return "isC_MU";
                case Type::isC_E:
                    return "isC_E";
                case Type::isC:
                    return "isC";
                case Type::isCC:
                    return "isCC";
                case Type::isGCC:
                    return "isGCC";
                case Type::isS:
                    return "isS"; 
                case Type::isS_MU:
                    return "isS_MU";
                case Type::isS_E:
                    return "isS_E";               
                case Type::isUD:
                    return "isUD";
                case Type::isUD_MU:
                    return "isUD_MU";
                case Type::isUD_E:
                    return "isUD_E";
                case Type::isG:
                    return "isG";
                case Type::isG_MU:
                    return "isG_MU";
                case Type::isG_E:
                    return "isG_E";
                
                case Type::isLLP_RAD:
                    return "isLLP_RAD";
                    
                case Type::isLLP_MU:
                    return "isLLP_MU";
                case Type::isLLP_E:
                    return "isLLP_E";
                    
                case Type::isLLP_Q:
                    return "isLLP_Q";
                    
                case Type::isLLP_QMU:
                    return "isLLP_QMU";
                case Type::isLLP_QE:
                    return "isLLP_QE";
                    
                case Type::isLLP_QQ:
                    return "isLLP_QQ";
                    
                case Type::isLLP_QQMU:
                    return "isLLP_QQMU";
                case Type::isLLP_QQE:
                    return "isLLP_QQE"; 
                    
                case Type::isLLP_B:
                    return "isLLP_B";
                    
                case Type::isLLP_BMU:
                    return "isLLP_BMU";
                case Type::isLLP_BE:
                    return "isLLP_BE";
                    
                case Type::isLLP_BB:
                    return "isLLP_BB";
                    
                case Type::isLLP_BBMU:
                    return "isLLP_BBMU";
                case Type::isLLP_BBE:
                    return "isLLP_BBE";

                case Type::isLLP_TAU:
                    return "isLLP_TAU";
                case Type::isLLP_QTAU:
                    return "isLLP_QTAU";
                case Type::isLLP_QQTAU:
                    return "isLLP_QQTAU";
                    
                case Type::isUndefined:
                    return "isUndefined";
            }
            return "isUndefined";
        }
};



}

#endif
