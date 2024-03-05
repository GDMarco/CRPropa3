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
    
    std::vector<double> getabEnergy() { //Isotropic
        return tabEnergy;
    }
    std::vector<double> getabRate() {
        return tabRate;
    }
    std::vector<double> getabE() {
        return tabE;
    }
    std::vector<double> getabs() {
        return tabs;
    }
    std::vector<std::vector<double>> getabCDF() {
        return tabCDF;
    }
    
    void setabEnergy (std::vector<double>& newtabEnergy) {
        tabEnergy = newtabEnergy;
    }
    void setabRate (std::vector<double>& newtabRate) {
        tabRate = newtabRate;
    }
    void setabE (std::vector<double>& newtabE) {
        tabE = newtabE;
    }
    void setabs (std::vector<double>& newtabs) {
        tabs = newtabs;
    }
    void setabCDF (std::vector<std::vector<double>>& newtabCDF) {
        tabCDF = newtabCDF;
    }
    
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
    
    std::vector<std::vector<double>> getabEnergy() { //PositionDependent
        return tabEnergy;
    }
    std::vector<std::vector<double>> getabRate() {
        return tabRate;
    }
    std::vector<std::vector<double>> getabE() {
        return tabE;
    }
    std::vector<std::vector<double>> getabs() {
        return tabs;
    }
    std::vector<std::vector<std::vector<double>>> getabCDF() {
        return tabCDF;
    }
    std::unordered_map<int, Vector3d> getphotonDict() {
        return photonDict;
    }
    
    void setabEnergy (std::vector<std::vector<double>>& newtabEnergy) {
        tabEnergy = newtabEnergy;
    }
    void setabRate (std::vector<std::vector<double>>& newtabRate) {
        tabRate = newtabRate;
    }
    void setabE (std::vector<std::vector<double>>& newtabE) {
        tabE = newtabE;
    }
    void setabs (std::vector<std::vector<double>>& newtabs) {
        tabs = newtabs;
    }
    void setabCDF (std::vector<std::vector<std::vector<double>>>& newtabCDF) {
        tabCDF = newtabCDF;
    }
    void setphotonDict (std::unordered_map<int, Vector3d>& newphotonDict) {
        photonDict = newphotonDict;
    }


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

//default constructors
InteractionRatesIsotropic::InteractionRatesIsotropic() : InteractionRates() { };

InteractionRatesPositionDependent::InteractionRatesPositionDependent() : InteractionRates() { };

// maybe I need to initialise the functions in each!

} // namespace crpropa

#endif // CRPROPA_INTERACTIONRATES_H
