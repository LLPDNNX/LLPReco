#include "LLPReco/DataFormats/interface/DisplacedGenVertex.h"

namespace llpdnnx
{

double DisplacedGenVertex::d3d() const
{
    return std::sqrt((hardInteraction-vertex).mag2());
}

double DisplacedGenVertex::dx() const
{
    return std::fabs(hardInteraction.x()-vertex.x());
}

double DisplacedGenVertex::dy() const
{
    return std::fabs(hardInteraction.y()-vertex.y());
}

double DisplacedGenVertex::dz() const
{
    return std::fabs(hardInteraction.z()-vertex.z());
}

double DisplacedGenVertex::dxy() const
{
    const double x = dx();
    const double y = dy();
    return std::sqrt(x*x+y*y);
}

}
