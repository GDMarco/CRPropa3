#ifndef CRPROPA_INTERACTIONRATES_H
#define CRPROPA_INTERACTIONRATES_H

#include "crpropa/Common.h"
#include "crpropa/Referenced.h"
#include "crpropa/Vector3.h"

#include <vector>
#include <string>
#include <unordered_map>

namespace crpropa {
/**
 * \addtogroup InteractionRates
 * @{
 */

/**
 @class Interaction Rates
 @brief Abstract base class for photon fields.
 */
class InteractionRates: public Referenced {
public:
    InteractionRates() {
}
    
};

/**
 @class InteractionRateIsotropic
 @brief Interaction rates decorator for tabulated isotropic interaction rates.
 */
class InteractionRatesIsotropic: public InteractionRates {
public:
    
    InteractionRatesIsotropic();
    
    std::vector<double> getabEnergy();
    std::vector<double> getabRate();
    std::vector<double> getabE();
    std::vector<double> getabs();
    std::vector<std::vector<double>> getabCDF();
    
    void setabEnergy (std::vector<double>& newtabEnergy);
    void setabRate (std::vector<double>& newtabRate);
    void setabE (std::vector<double>& newtabE);
    void setabs (std::vector<double>& newtabs);
    void setabCDF (std::vector<std::vector<double>>& newtabCDF);
    
protected:
    
    // tabulated interaction rates 1/lambda(E)
    std::vector<double> tabEnergy; //!< electron energy in [J]
    std::vector<double> tabRate; //!< interaction rate in [1/m]
    
    // tabulated CDF(s_kin, E) = cumulative differential interaction rate
    std::vector<double> tabE; //!< electron energy in [J]
    std::vector<double> tabs; //!< s_kin = s - m^2 in [J**2]
    std::vector<std::vector<double>> tabCDF; //!< cumulative interaction rate
    
};

class InteractionRatesPositionDependent: public InteractionRates {
public:
    
    InteractionRatesPositionDependent();
    
    std::vector<std::vector<double>> getabEnergy();
    std::vector<std::vector<double>> getabRate();
    std::vector<std::vector<double>> getabE();
    std::vector<std::vector<double>> getabs();
    std::vector<std::vector<std::vector<double>>> getabCDF();
    std::unordered_map<int, Vector3d> getphotonDict();
    
    void setabEnergy (std::vector<std::vector<double>>& newtabEnergy);
    void setabRate (std::vector<std::vector<double>>& newtabRate);
    void setabE (std::vector<std::vector<double>>& newtabE);
    void setabs (std::vector<std::vector<double>>& newtabs);
    void setabCDF (std::vector<std::vector<std::vector<double>>>& newtabCDF);
    void setphotonDict (std::unordered_map<int, Vector3d>& newphotonDict);


protected:
    
    // tabulated interaction rates 1/lambda(E)
    std::vector<std::vector<double>> tabEnergy;
    std::vector<std::vector<double>> tabRate;
    
    // tabulated CDF(s_kin, E) = cumulative differential interaction rate
    std::vector<std::vector<double>> tabE;
    std::vector<std::vector<double>> tabs;
    std::vector<std::vector<std::vector<double>>> tabCDF;
    std::unordered_map<int, Vector3d> photonDict;
};

} // namespace crpropa

#endif // CRPROPA_INTERACTIONRATES_H
