#ifndef NANOX_LLPSPECIFIC_LLPINFO_H
#define NANOX_LLPSPECIFIC_LLPINFO_H

#include "NANOX/DataFormats/interface/TagData.h"
#include <iostream>

namespace nanox
{

class LLPInfo:
    public TagData
{
    public:
        
    
    
        class Data:
            public PropertyContainer
        {
            public:
                float llp_mass;
                float llp_pt;
                float llp_eta;
                float llp_phi;
                
                float lsp_mass;
                float lsp_pt;
                float lsp_eta;
                float lsp_phi;
                
                int quarkFlavor;

                Data():
                    llp_mass(0),
                    llp_pt(0),
                    llp_eta(0),
                    llp_phi(0),
                    
                    lsp_mass(0),
                    lsp_pt(0),
                    lsp_eta(0),
                    lsp_phi(0),
                    
                    quarkFlavor(0)
                {
                }
                
        };
        std::vector<Data> llpData;

        virtual void saveTagData(ArchiveInterface& archive) const override
        {
            ArrayInterface& dataArray = archive.initArray("llpinfo",llpData.size());
            
            dataArray.bookProperty("llp_mass",&Data::llp_mass);
            dataArray.bookProperty("llp_pt",&Data::llp_pt);
            dataArray.bookProperty("llp_eta",&Data::llp_eta);
            dataArray.bookProperty("llp_phi",&Data::llp_phi);
            
            dataArray.bookProperty("lsp_mass",&Data::lsp_mass);
            dataArray.bookProperty("lsp_pt",&Data::lsp_pt);
            dataArray.bookProperty("lsp_eta",&Data::lsp_eta);
            dataArray.bookProperty("lsp_phi",&Data::lsp_phi);
            
            dataArray.bookProperty("quarkFlavor",&Data::quarkFlavor);
            
            for (unsigned int illp = 0; illp < llpData.size(); ++illp)
            {
                dataArray.fill(&llpData[illp],illp);
            }
        }
        
        virtual ~LLPInfo()
        {
        }
};

}

#endif
