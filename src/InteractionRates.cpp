#include "crpropa/InteractionRates.h"
#include "crpropa/Referenced.h"
#include "crpropa/Vector3.h"

#include <vector>
#include <string>
#include <unordered_map>

namespace crpropa {

InteractionRatesIsotropic::InteractionRatesIsotropic() : InteractionRates() { }

std::vector<double> InteractionRatesIsotropic::getabEnergy() {
    return tabEnergy;
}

std::vector<double> InteractionRatesIsotropic::getabRate() {
    return tabRate;
}

std::vector<double> InteractionRatesIsotropic::getabE() {
    return tabE;
}
std::vector<double> InteractionRatesIsotropic::getabs() {
    return tabs;
}
std::vector<std::vector<double>> InteractionRatesIsotropic::getabCDF() {
    return tabCDF;
}

void InteractionRatesIsotropic::setabEnergy (std::vector<double>& newtabEnergy) {
    tabEnergy = newtabEnergy;
}
void InteractionRatesIsotropic::setabRate (std::vector<double>& newtabRate) {
    tabRate = newtabRate;
}
void InteractionRatesIsotropic::setabE (std::vector<double>& newtabE) {
    tabE = newtabE;
}
void InteractionRatesIsotropic::setabs (std::vector<double>& newtabs) {
    tabs = newtabs;
}
void InteractionRatesIsotropic::setabCDF (std::vector<std::vector<double>>& newtabCDF) {
    tabCDF = newtabCDF;
}

InteractionRatesPositionDependent::InteractionRatesPositionDependent() : InteractionRates() {}

std::vector<std::vector<double>> InteractionRatesPositionDependent::getabEnergy() {
    return tabEnergy;
}
std::vector<std::vector<double>> InteractionRatesPositionDependent::getabRate() {
    return tabRate;
}
std::vector<std::vector<double>> InteractionRatesPositionDependent::getabE() {
    return tabE;
}
std::vector<std::vector<double>> InteractionRatesPositionDependent::getabs() {
    return tabs;
}
std::vector<std::vector<std::vector<double>>> InteractionRatesPositionDependent::getabCDF() {
    return tabCDF;
}
std::unordered_map<int, Vector3d> InteractionRatesPositionDependent::getphotonDict() {
    return photonDict;
}

void InteractionRatesPositionDependent::setabEnergy (std::vector<std::vector<double>>& newtabEnergy) {
    tabEnergy = newtabEnergy;
}
void InteractionRatesPositionDependent::setabRate (std::vector<std::vector<double>>& newtabRate) {
    tabRate = newtabRate;
}
void InteractionRatesPositionDependent::setabE (std::vector<std::vector<double>>& newtabE) {
    tabE = newtabE;
}
void InteractionRatesPositionDependent::setabs (std::vector<std::vector<double>>& newtabs) {
    tabs = newtabs;
}
void InteractionRatesPositionDependent::setabCDF (std::vector<std::vector<std::vector<double>>>& newtabCDF) {
    tabCDF = newtabCDF;
}
void InteractionRatesPositionDependent::setphotonDict (std::unordered_map<int, Vector3d>& newphotonDict) {
    photonDict = newphotonDict;
}

} //namespace crpropa
