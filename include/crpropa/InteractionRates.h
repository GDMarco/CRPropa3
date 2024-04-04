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
 @brief Abstract base class for photon fields interaction rates.
 */
class InteractionRates: public Referenced {
public:
    InteractionRates() {
        this->ratesName = "AbstractInteractionRates";
        this->isPositionDependent = false;
    }
    
    std::string getRatesName() const {
        return this->ratesName;
    }
    
    bool hasPositionDependence() const {
        return this->isPositionDependent;
    }
    
    void setRatesName(std::string ratesName) {
        this->ratesName = ratesName;
    }

protected: 

  std::string ratesName;
  bool isPositionDependent; 

};

/**
 @class InteractionRateHomogeneous
 @brief Interaction rates decorator for tabulated homogeneous interaction rates.
 */
class InteractionRatesHomogeneous: public InteractionRates {
public:
    InteractionRatesHomogeneous(const std::string ratesName, const bool isPositionDependent = true);
    
    std::vector<double> getTabulatedEnergy();
    std::vector<double> getTabulatedRate();
    std::vector<double> getTabulatedE();
    std::vector<double> getTabulateds();
    std::vector<std::vector<double>> getTabulatedCDF();
    
    void setTabulatedEnergy (std::vector<double>& tabEnergy);
    void setTabulatedRate (std::vector<double>& tabRate);
    void setTabulatedE (std::vector<double>& tabE);
    void setTabulateds (std::vector<double>& tabs);
    void setTabulatedCDF (std::vector<std::vector<double>>& tabCDF);
    
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
    InteractionRatesPositionDependent(const std::string ratesName, const bool isPositionDependent = true);
    
    std::vector<std::vector<double>> getTabulatedEnergy();
    std::vector<std::vector<double>> getTabulatedRate();
    std::vector<std::vector<double>> getTabulatedE();
    std::vector<std::vector<double>> getTabulateds();
    std::vector<std::vector<std::vector<double>>> getTabulatedCDF();
    std::unordered_map<int, Vector3d> getPhotonDict();
    
    void setTabulatedEnergy (std::vector<std::vector<double>>& tabEnergy);
    void setTabulatedRate (std::vector<std::vector<double>>& tabRate);
    void setTabulatedE (std::vector<std::vector<double>>& tabE);
    void setTabulateds (std::vector<std::vector<double>>& tabs);
    void setTabulatedCDF (std::vector<std::vector<std::vector<double>>>& tabCDF);
    void setPhotonDict (std::unordered_map<int, Vector3d>& photonDict); 

protected:
    
    // tabulated interaction rates 1/lambda(E)
    std::vector<std::vector<double>> tabEnergy; //!< electron energy in [J]
    std::vector<std::vector<double>> tabRate; //!< interaction rate in [1/m]
    
    // tabulated CDF(s_kin, E) = cumulative differential interaction rate
    std::vector<std::vector<double>> tabE; //!< electron energy in [J]
    std::vector<std::vector<double>> tabs; //!< s_kin = s - m^2 in [J**2]
    std::vector<std::vector<std::vector<double>>> tabCDF; //!< cumulative interaction rate
    std::unordered_map<int, Vector3d> photonDict; //!< dictionary to link tables to spatial coordinates
};

} // namespace crpropa

#endif // CRPROPA_INTERACTIONRATES_H
