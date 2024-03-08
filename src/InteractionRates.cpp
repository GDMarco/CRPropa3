#include "crpropa/InteractionRates.h"
#include "crpropa/Referenced.h"
#include "crpropa/Vector3.h"

#include <vector>
#include <string>
#include <unordered_map>

namespace crpropa {

InteractionRatesIsotropic::InteractionRatesIsotropic(std::string ratesName, bool isPositionDependent) : InteractionRates() {
  this->ratesName = ratesName;
  this->isPositionDependent = isPositionDependent;
}

std::vector<double> InteractionRatesIsotropic::getTabulatedEnergy() {
    return tabEnergy;
}

std::vector<double> InteractionRatesIsotropic::getTabulatedRate() {
    return tabRate;
}

std::vector<double> InteractionRatesIsotropic::getTabulatedE() {
    return tabE;
}
std::vector<double> InteractionRatesIsotropic::getTabulateds() {
    return tabs;
}
std::vector<std::vector<double>> InteractionRatesIsotropic::getTabulatedCDF() {
    return tabCDF;
}

void InteractionRatesIsotropic::setTabulatedEnergy (std::vector<double>& tabEnergy) {
    tabEnergy = tabEnergy;
}
void InteractionRatesIsotropic::setTabulatedRate (std::vector<double>& tabRate) {
    tabRate = tabRate;
}
void InteractionRatesIsotropic::setTabulatedE (std::vector<double>& tabE) {
    tabE = tabE;
}
void InteractionRatesIsotropic::setTabulateds (std::vector<double>& tabs) {
    tabs = tabs;
}
void InteractionRatesIsotropic::setTabulatedCDF (std::vector<std::vector<double>>& tabCDF) {
    tabCDF = tabCDF;
}

InteractionRatesPositionDependent::InteractionRatesPositionDependent(std::string ratesName, bool isPositionDependent) : InteractionRates() {
  this->ratesName = ratesName;
  this->isPositionDependent = isPositionDependent;
}

std::vector<std::vector<double>> InteractionRatesPositionDependent::getTabulatedEnergy() {
    return tabEnergy;
}
std::vector<std::vector<double>> InteractionRatesPositionDependent::getTabulatedRate() {
    return tabRate;
}
std::vector<std::vector<double>> InteractionRatesPositionDependent::getTabulatedE() {
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

void InteractionRatesPositionDependent::setTabulatedEnergy (std::vector<std::vector<double>>& tabEnergy) {
    tabEnergy = tabEnergy;
}
void InteractionRatesPositionDependent::setTabulatedRate (std::vector<std::vector<double>>& tabRate) {
    tabRate = tabRate;
}
void InteractionRatesPositionDependent::setTabulatedE (std::vector<std::vector<double>>& tabE) {
    tabE = tabE;
}
void InteractionRatesPositionDependent::setTabulateds (std::vector<std::vector<double>>& tabs) {
    tabs = tabs;
}
void InteractionRatesPositionDependent::setTabulatedCDF (std::vector<std::vector<std::vector<double>>>& tabCDF) {
    tabCDF = tabCDF;
}
void InteractionRatesPositionDependent::setPhotonDict (std::unordered_map<int, Vector3d>& photonDict) {
    photonDict = photonDict;
}

} //namespace crpropa
