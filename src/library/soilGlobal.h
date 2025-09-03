#ifndef SOILGLOBAL_H
#define SOILGLOBAL_H

#include "fvCFD.H"

#ifdef VSSFOAM_DEBUG
#define DEBUG(x)         \
    do                   \
    {                    \
        Info << x << nl; \
    } while (0)
#else
#define DEBUG(x) \
    do           \
    {            \
    } while (0)
#endif

namespace Soil::Global
{

    void checkUnits(const dimensionedScalar &arg, const dimensionSet &expectedUnits, const std::string &argName);
    void checkUnits(const volScalarField &arg, const dimensionSet &expectedUnits, const std::string &argName);

}

#endif // SOILGLOBAL_H
