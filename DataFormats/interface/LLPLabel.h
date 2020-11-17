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
            isLeptonic_B,
            isC,
            isCC,
            isLeptonic_C,
            isS,
            isUD,
            isG,
            isPrompt_MU,
            isPrompt_E,
            isPrompt_PHOTON,
            isPrompt_TAU,
            
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
            isLLP_BTAU,
            isLLP_BBTAU,
            
            isLLP_PHOTON,
            isLLP_QPHOTON,
            isLLP_QQPHOTON,
            isLLP_BPHOTON,
            isLLP_BBPHOTON,
            isUndefined
            
        };
        
        Type type;
        
        enum class TauDecay {
            NO_TAU,     //no tau decay
            INVISIBLE,  //tau decay but not reconstructable
            E,          //to electron
            MU,         //to muon
            H,          //1 charged hadron
            H_1PI0,     //1 charged hadron + pi0(->2gamma)
            H_XPI0,     //1 charged hadron + 2 or more pi0(->2gamma)
            HHH,         //3 charged hadrons
            HHH_XPI0     //3 charged hadron + 1 or more pi0(->2gamma)
        };
        
        TauDecay tauDecay;
        
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
        float sharedVertexFraction;
        
        float genTauMass;
        float recoTauMass;
        
        LLPLabel():
            type(Type::isUndefined),
            tauDecay(TauDecay::NO_TAU),
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
            matchedGenJetPt(-1),
            sharedVertexFraction(0),
            genTauMass(-1),
            recoTauMass(-1)
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
                case Type::isLeptonic_B:
                    return "isLeptonic_B";
                case Type::isLeptonic_C:
                    return "isLeptonic_C";
                case Type::isC:
                    return "isC";
                case Type::isCC:
                    return "isCC";
                case Type::isS:
                    return "isS";               
                case Type::isUD:
                    return "isUD";
                case Type::isG:
                    return "isG";

                case Type::isPrompt_MU:
                    return "isPrompt_MU";
                case Type::isPrompt_E:
                    return "isPrompt_E";
                case Type::isPrompt_TAU:
                    return "isPrompt_TAU";
                case Type::isPrompt_PHOTON:
                    return "isPrompt_PHOTON";
                
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
                case Type::isLLP_BTAU:
                    return "isLLP_BTAU";
                case Type::isLLP_BBTAU:
                    return "isLLP_BBTAU";
                    
                case Type::isLLP_PHOTON:
                    return "isLLP_TAU";
                case Type::isLLP_QPHOTON:
                    return "isLLP_QTAU";
                case Type::isLLP_QQPHOTON:
                    return "isLLP_QQTAU";
                case Type::isLLP_BPHOTON:
                    return "isLLP_BTAU";
                case Type::isLLP_BBPHOTON:
                    return "isLLP_BBTAU";
                    
                case Type::isUndefined:
                    return "isUndefined";
            }
            return "isUndefined";
        }
};



}

#endif
