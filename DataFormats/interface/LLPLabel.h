#ifndef LLPReco_DataFormats_LLPLabel_h
#define LLPReco_DataFormats_LLPLabel_h

#include "DataFormats/BTauReco/interface/FeaturesTagInfo.h"

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
            isC,
            isCC,
            isGCC,
            isS,
            isUD,
            isG,
            
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
            isUndefined
            
        };
        
        Type type;
        
        int jetIdx;
        int partonFlavor;
        int hadronFlavor;
        int llpId;

        float displacement;
        float displacement_xy;
        float displacement_z;
        float decay_angle;
        float betagamma;
        
        
        LLPLabel():
            type(Type::isUndefined),
            partonFlavor(0),
            hadronFlavor(0),
            llpId(0),
            displacement(-10), //log10(x/1cm)
            displacement_xy(-10), //log10(x/1cm)
            displacement_z(-10), //log10(x/1cm)
            decay_angle(0),
            betagamma(0)
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
                case Type::isC:
                    return "isC";
                case Type::isCC:
                    return "isCC";
                case Type::isGCC:
                    return "isGCC";
                case Type::isS:
                    return "isS";                
                case Type::isUD:
                    return "isUD";
                case Type::isG:
                    return "isG";
                
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
                    
                case Type::isUndefined:
                    return "isUndefined";
            }
            return "isUndefined";
        }
};



}

#endif
