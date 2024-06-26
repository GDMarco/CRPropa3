#include "crpropa/InteractionRates.h"
#include "crpropa/Referenced.h"
#include "crpropa/Vector3.h"
#include "crpropa/Units.h"

#include <vector>
#include <string>
#include <unordered_map>

namespace crpropa {

InteractionRatesHomogeneous::InteractionRatesHomogeneous(std::string ratesName, bool isPositionDependent) : InteractionRates() {
  this->ratesName = ratesName;
  this->isPositionDependent = isPositionDependent;
}

std::vector<double> InteractionRatesHomogeneous::getTabulatedEnergy() {
    return tabEnergy;
}

std::vector<double> InteractionRatesHomogeneous::getTabulatedRate() {
    return tabRate;
}

std::vector<double> InteractionRatesHomogeneous::getTabulatedE() {
    return tabE;
}
std::vector<double> InteractionRatesHomogeneous::getTabulateds() {
    return tabs;
}
std::vector<std::vector<double>> InteractionRatesHomogeneous::getTabulatedCDF() {
    return tabCDF;
}

void InteractionRatesHomogeneous::setTabulatedEnergy (std::vector<double>& tabEnergy) {
    this->tabEnergy = tabEnergy;
}
void InteractionRatesHomogeneous::setTabulatedRate (std::vector<double>& tabRate) {
    this->tabRate = tabRate;
}
void InteractionRatesHomogeneous::setTabulatedE (std::vector<double>& tabE) {
    this->tabE = tabE;
}
void InteractionRatesHomogeneous::setTabulateds (std::vector<double>& tabs) {
    this->tabs = tabs;
}
void InteractionRatesHomogeneous::setTabulatedCDF (std::vector<std::vector<double>>& tabCDF) {
    this->tabCDF = tabCDF;
}

InteractionRatesPositionDependent::InteractionRatesPositionDependent(std::string ratesName, bool isPositionDependent) : InteractionRates() {
  this->ratesName = ratesName;
  this->isPositionDependent = isPositionDependent;
}

std::vector<double> InteractionRatesPositionDependent::getTabulatedEnergy() {
    return tabEnergy;
}
std::vector<std::vector<double>> InteractionRatesPositionDependent::getTabulatedRate() {
    return tabRate;
}
std::vector<double> InteractionRatesPositionDependent::getTabulatedE() {
    return tabE;
}
std::vector<std::vector<double>> InteractionRatesPositionDependent::getTabulateds() {
    return tabs;
}
std::vector<std::vector<std::vector<double>>> InteractionRatesPositionDependent::getTabulatedCDF() {
    return tabCDF;
}
std::unordered_map<int, Vector3d> InteractionRatesPositionDependent::getPhotonDict() {
    return photonDict;
}

std::vector<double> InteractionRatesPositionDependent::getClosestRate(const Vector3d &position) {
    
    std::unordered_map<int,Vector3d> photonDict = this->getPhotonDict();
    
    double dMin = 1000. * kpc;
    int iMin = -1;
    
    for (const auto& el : photonDict) {
        
        Vector3d posNode = el.second;
        double d;
        d = sqrt((- posNode.x / kpc - position.x / kpc) * (- posNode.x / kpc - position.x / kpc) + (posNode.y / kpc - position.y / kpc) * (posNode.y / kpc - position.y / kpc) + (posNode.z / kpc - position.z / kpc) * (posNode.z / kpc - position.z / kpc));
        
        if (d < dMin) {
            dMin = d;
            iMin = el.first;
        }
    }
    return tabRate[iMin];
}

std::vector<double> InteractionRatesPositionDependent::getClosests(const Vector3d &position) {
    
    std::unordered_map<int,Vector3d> photonDict = this->getPhotonDict();
    
    double dMin = 1000. * kpc;
    int iMin = -1;
    
    for (const auto& el : photonDict) {
        
        Vector3d posNode = el.second;
        double d;
        d = sqrt((- posNode.x / kpc - position.x / kpc) * (- posNode.x / kpc - position.x / kpc) + (posNode.y / kpc - position.y / kpc) * (posNode.y / kpc - position.y / kpc) + (posNode.z / kpc - position.z / kpc) * (posNode.z / kpc - position.z / kpc));
        
        if (d < dMin) {
            dMin = d;
            iMin = el.first;
        }
    }
    
    return tabs[iMin];
}

std::vector<std::vector<double>> InteractionRatesPositionDependent::getClosestCDF(const Vector3d &position) {
    
    std::unordered_map<int,Vector3d> photonDict = this->getPhotonDict();
    
    double dMin = 1000. * kpc;
    int iMin = -1;
    
    for (const auto& el : photonDict) {
        
        Vector3d posNode = el.second;
        double d;
        d = sqrt((- posNode.x / kpc - position.x / kpc) * (- posNode.x / kpc - position.x / kpc) + (posNode.y / kpc - position.y / kpc) * (posNode.y / kpc - position.y / kpc) + (posNode.z / kpc - position.z / kpc) * (posNode.z / kpc - position.z / kpc));
        
        if (d < dMin) {
            dMin = d;
            iMin = el.first;
        }
    }
    return tabCDF[iMin];
}

void InteractionRatesPositionDependent::setTabulatedEnergy (std::vector<double>& tabEnergy) {
    this->tabEnergy = tabEnergy;
}
void InteractionRatesPositionDependent::setTabulatedRate (std::vector<std::vector<double>>& tabRate) {
    this->tabRate = tabRate;
}
void InteractionRatesPositionDependent::setTabulatedE (std::vector<double>& tabE) {
    this->tabE = tabE;
}
void InteractionRatesPositionDependent::setTabulateds (std::vector<std::vector<double>>& tabs) {
    this->tabs = tabs;
}
void InteractionRatesPositionDependent::setTabulatedCDF (std::vector<std::vector<std::vector<double>>>& tabCDF) {
    this->tabCDF = tabCDF;
}
void InteractionRatesPositionDependent::setPhotonDict (std::unordered_map<int, Vector3d>& photonDict) {
    this->photonDict = photonDict;
}

} //namespace crpropa
