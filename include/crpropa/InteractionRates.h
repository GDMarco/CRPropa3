#ifndef CRPROPA_INTERACTIONRATES_H
#define CRPROPA_INTERACTIONRATES_H

#include "crpropa/Common.h"
#include "crpropa/Referenced.h"
#include "crpropa/Vector3.h"

#include <nanoflann.hpp>

#include <vector>
#include <string>
#include <unordered_map>

namespace crpropa {

struct PointCloud {
    std::vector<Vector3d> points;
    std::vector<int> ids;

    inline size_t kdtree_get_point_count() const { return points.size(); }

    // Returns the X/Y/Z coordinate of the ith point
    inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
        if (dim == 0) return points[idx].x;
        if (dim == 1) return points[idx].y;
        return points[idx].z;
    }

    // Optional bounding-box computation (required by nanoflann)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const {
        return false;  // No bounding box optimization
    }
    
};

// Define KD-Tree index type
using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, PointCloud>,
    PointCloud,
    3  // 3D space
>;


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
    
    virtual double getProcessRate(const double E, const Vector3d &position) const = 0;
    virtual void loadPerformInteractionTabs(const Vector3d &position, std::vector<double> &tabE, std::vector<double> &tabs, std::vector<std::vector<double>> &tabCDF) const = 0;
    
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
    
    std::vector<double> getTabulatedEnergy() const;
    std::vector<double> getTabulatedRate() const;
    std::vector<double> getTabulatedE() const;
    std::vector<double> getTabulateds() const;
    std::vector<std::vector<double>> getTabulatedCDF() const;
    
    double getProcessRate(const double E, const Vector3d &position) const;
    void loadPerformInteractionTabs(const Vector3d &position, std::vector<double> &tabE, std::vector<double> &tabs, std::vector<std::vector<double>> &tabCDF) const;
    
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
    
    int findClosestGridPoint(const Vector3d &position) const; 
    
    std::vector<double> getTabulatedEnergy() const;
    std::vector<std::vector<double>> getTabulatedRate() const;
    std::vector<double> getTabulatedE() const;
    std::vector<std::vector<double>> getTabulateds() const;
    std::vector<std::vector<std::vector<double>>> getTabulatedCDF() const;
    std::unordered_map<int, Vector3d> getPhotonDict() const;
    std::vector<double> getClosestRate(const Vector3d &position) const;
    std::vector<double> getClosests(const Vector3d &position) const;
    std::vector<std::vector<double>> getClosestCDF(const Vector3d &position) const;
    
    double getProcessRate(const double E, const Vector3d &position) const;
    void loadPerformInteractionTabs(const Vector3d &position, std::vector<double> &tabE, std::vector<double> &tabs, std::vector<std::vector<double>> &tabCDF) const;
    
    void setTabulatedEnergy (std::vector<double>& tabEnergy);
    void setTabulatedRate (std::vector<std::vector<double>>& tabRate);
    void setTabulatedE (std::vector<double>& tabE);
    void setTabulateds (std::vector<std::vector<double>>& tabs);
    void setTabulatedCDF (std::vector<std::vector<std::vector<double>>>& tabCDF);
    void setPhotonDict (std::unordered_map<int, Vector3d>& photonDict);

protected:
    
    // tabulated interaction rates 1/lambda(E)
    std::vector<double> tabEnergy; //!< electron energy in [J], assuming the same energy binning in each node
    std::vector<std::vector<double>> tabRate; //!< interaction rate in [1/m]
    
    // tabulated CDF(s_kin, E) = cumulative differential interaction rate
    std::vector<double> tabE; //!< electron energy in [J], assuming the same energy binning in each node
    std::vector<std::vector<double>> tabs; //!< s_kin = s - m^2 in [J**2]
    std::vector<std::vector<std::vector<double>>> tabCDF; //!< cumulative interaction rate
    std::unordered_map<int, Vector3d> photonDict; //!< dictionary to link tables to spatial coordinates
    
    
    PointCloud cloud; //!< Point cloud for nanoflann KD tree
    KDTree* tree = nullptr; //!< Pointer to the KD-Tree//!
    
};

} // namespace crpropa

#endif // CRPROPA_INTERACTIONRATES_H
