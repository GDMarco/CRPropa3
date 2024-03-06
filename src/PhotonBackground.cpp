#include "crpropa/PhotonBackground.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/Vector3.h"

#include "kiss/logger.h"

#include <fstream>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <filesystem>
#include <string>
#include <sstream>
#include <unordered_map>

namespace crpropa {

TabularPhotonField::TabularPhotonField(std::string fieldName, bool isRedshiftDependent, bool isSpatialDependent) {
	this->fieldName = fieldName;
	this->isRedshiftDependent = isRedshiftDependent;
    this->isSpatialDependent = isSpatialDependent;

    if (this->isSpatialDependent) {
        KISS_LOG_WARNING << "Photon Field " << fieldName << " is spatial dependent! It is not the correct class. \n";
        exit(1); // to check! 
    } else {
        readPhotonEnergy(getDataPath("") + "Scaling/" + this->fieldName + "_photonEnergy.txt");
        readPhotonDensity(getDataPath("") + "Scaling/" + this->fieldName + "_photonDensity.txt");
        if (this->isRedshiftDependent)
            readRedshift(getDataPath("") + "Scaling/" + this->fieldName + "_redshift.txt");
        
        checkInputData();
        
        if (this->isRedshiftDependent)
            initRedshiftScaling();
    }
}


double TabularPhotonField::getPhotonDensity(double Ephoton, double z, const Vector3d &pos) const {
	if (this->isRedshiftDependent) {
		// fix behaviour for future redshift. See issue #414
		// with redshift < 0 the photon density is set to 0 in interpolate2d. 
		// Therefore it is assumed that the photon density does not change from values at z = 0. This is only valid for small changes in redshift.
		double zMin = this->redshifts[0];
		if(z < zMin){
			if(z < -1) {
				KISS_LOG_WARNING << "Photon Field " << fieldName << " uses FutureRedshift with z < -1. The photon density is set to n(Ephoton, z=0). \n";
			}
			return getPhotonDensity(Ephoton, zMin);
		} else {
			return interpolate2d(Ephoton, z, this->photonEnergies, this->redshifts, this->photonDensity);
		}
	} else {
		return interpolate(Ephoton, this->photonEnergies, this->photonDensity);
	}
}

double TabularPhotonField::getRedshiftScaling(double z) const {
	if (!this->isRedshiftDependent)
		return 1.;
 
	if (z < this->redshifts.front())
		return 1.;
 
	if (z > this->redshifts.back())
		return 0.;
 
	return interpolate(z, this->redshifts, this->redshiftScalings);
}

double TabularPhotonField::getMinimumPhotonEnergy(double z, const Vector3d &pos) const{
	return photonEnergies[0];
}

double TabularPhotonField::getMaximumPhotonEnergy(double z, const Vector3d &pos) const{
	return photonEnergies[photonEnergies.size() -1];
}

void TabularPhotonField::readPhotonEnergy(std::string filePath) {
	std::ifstream infile(filePath.c_str());
	if (!infile.good())
		throw std::runtime_error("TabularPhotonField::readPhotonEnergy: could not open " + filePath);

	std::string line;
	while (std::getline(infile, line)) {
		if ((line.size() > 0) & (line[0] != '#') )
			this->photonEnergies.push_back(std::stod(line));
	}
	infile.close();
}

void TabularPhotonField::readPhotonDensity(std::string filePath) {
	std::ifstream infile(filePath.c_str());
	if (!infile.good())
		throw std::runtime_error("TabularPhotonField::readPhotonDensity: could not open " + filePath);

	std::string line;
	while (std::getline(infile, line)) {
		if ((line.size() > 0) & (line[0] != '#') )
			this->photonDensity.push_back(std::stod(line));
	}
	infile.close();
}

void TabularPhotonField::readRedshift(std::string filePath) {
	std::ifstream infile(filePath.c_str());
	if (!infile.good())
		throw std::runtime_error("TabularPhotonField::initRedshift: could not open " + filePath);

	std::string line;
	while (std::getline(infile, line)) {
		if ((line.size() > 0) & (line[0] != '#') )
			this->redshifts.push_back(std::stod(line));
	}
	infile.close();
}

void TabularPhotonField::initRedshiftScaling() {
	double n0 = 0.;
	for (int i = 0; i < this->redshifts.size(); ++i) {
		double z = this->redshifts[i];
		double n = 0.;
		for (int j = 0; j < this->photonEnergies.size()-1; ++j) {
			double e_j = this->photonEnergies[j];
			double e_j1 = this->photonEnergies[j+1];
			double deltaLogE = std::log10(e_j1) - std::log10(e_j);
			if (z == 0.)
				n0 += (getPhotonDensity(e_j, 0) + getPhotonDensity(e_j1, 0)) / 2. * deltaLogE;
			n += (getPhotonDensity(e_j, z) + getPhotonDensity(e_j1, z)) / 2. * deltaLogE;
		}
		this->redshiftScalings.push_back(n / n0);
	}
}

void TabularPhotonField::checkInputData() const {
	if (this->isRedshiftDependent) {
		if (this->photonDensity.size() != this->photonEnergies.size() * this-> redshifts.size())
			throw std::runtime_error("TabularPhotonField::checkInputData: length of photon density input is unequal to length of photon energy input times length of redshift input");
	} else {
		if (this->photonEnergies.size() != this->photonDensity.size())
			throw std::runtime_error("TabularPhotonField::checkInputData: length of photon energy input is unequal to length of photon density input");
	}

	for (int i = 0; i < this->photonEnergies.size(); ++i) {
		double ePrevious = 0.;
		double e = this->photonEnergies[i];
		if (e <= 0.)
			throw std::runtime_error("TabularPhotonField::checkInputData: a value in the photon energy input is not positive");
		if (e <= ePrevious)
			throw std::runtime_error("TabularPhotonField::checkInputData: photon energy values are not strictly increasing");
		ePrevious = e;
	}

	for (int i = 0; i < this->photonDensity.size(); ++i) {
		if (this->photonDensity[i] < 0.)
			throw std::runtime_error("TabularPhotonField::checkInputData: a value in the photon density input is negative");
	}

	if (this->isRedshiftDependent) {
		if (this->redshifts[0] != 0.)
			throw std::runtime_error("TabularPhotonField::checkInputData: redshift input must start with zero");

		for (int i = 0; i < this->redshifts.size(); ++i) {
			double zPrevious = -1.;
			double z = this->redshifts[i];
			if (z < 0.)
				throw std::runtime_error("TabularPhotonField::checkInputData: a value in the redshift input is negative");
			if (z <= zPrevious)
				throw std::runtime_error("TabularPhotonField::checkInputData: redshift values are not strictly increasing");
			zPrevious = z;
		}

		for (int i = 0; i < this->redshiftScalings.size(); ++i) {
			double scalingFactor = this->redshiftScalings[i];
			if (scalingFactor <= 0.)
				throw std::runtime_error("TabularPhotonField::checkInputData: initRedshiftScaling has created a non-positive scaling factor");
		}
	}
}

/**
Tabular spatial photon field, _photonEnergy.txt and _photonDensity.txt files have to be into different repositories.
*/
TabularSpatialPhotonField::TabularSpatialPhotonField(std::string fieldName, bool isRedshiftDependent, bool isSpatialDependent) {
    this->fieldName = fieldName;
    this->isRedshiftDependent = isRedshiftDependent;
    this->isSpatialDependent = isSpatialDependent;
    
    if (this->isRedshiftDependent) {
        
        KISS_LOG_WARNING << "Photon Field " << fieldName << " is redshift dependent! It is not the correct class. \n";
        exit(1); // to check
        
    } else if (!this->isSpatialDependent) { 
        
        KISS_LOG_WARNING << "Photon Field " << fieldName << " is not spatial dependent! It is not the correct class. \n";
        exit(1); // to check
        
    } else {
        std::cout << "enter TabPhotonField" << std::endl;
        std::__fs::filesystem::path dirE = getDataPath("") + "Scaling/" + this->fieldName + "/photonEnergy/";
        std::cout << "dirE TabPhotonField" << std::endl;
        std::unordered_map<int, Vector3d> photonDict;
        int iFile = 0;
        
        //for cycle for all the files in the path: build Dictionary + build this->photonEnergy
        for (auto const& dir_entry : std::__fs::filesystem::directory_iterator{dirE}) {
            
            double x, y, z;
            std::string str;
            std::stringstream ss;
            
            //Taking the filename
            std::string filename = splitFilename(dir_entry.path().string());
            ss << filename;
            
            //Getline function to take the x, y, z coordinates and computing the distance d
            int iLine = 0;
            
            while (getline(ss, str, '_')) {
                if (iLine == 2) {
                    x = stod(str) * kpc;
                }
                if (iLine == 3) {
                    y = stod(str) * kpc;
                }
                if (iLine == 4) {
                    z = stod(str) * kpc;
                }
                iLine = iLine + 1;
            }
            
            Vector3d vPos(x, y, z);
            photonDict[iFile] = vPos;
            
            iFile = iFile + 1;
            
            //charge all the fields in the proper double vector
            std::vector<double> vE = readPhotonEnergy(dir_entry.path().string());
            
            this->photonEnergies.push_back(vE);
            
        }
        
        this->photonDict = photonDict;
        
        std::__fs::filesystem::path dirD = getDataPath("") + "Scaling/" + this->fieldName + "/photonDensity/";
        
        for (auto const& dir_entry : std::__fs::filesystem::directory_iterator{dirD}) {
            
            std::vector<double> vD = readPhotonDensity(dir_entry.path().string());
            
            this->photonDensity.push_back(vD);
            
        }
        checkInputData();
    }
}

std::string TabularSpatialPhotonField::splitFilename(const std::string str) const {
    std::size_t found = str.find_last_of("/\\");
    std::string s = str.substr(found+1);
    return s;
}

double TabularSpatialPhotonField::getPhotonDensity(const double ePhoton, double z, const Vector3d &pos) const {
    
    double dMin = 1000.;
    int iMin = -1;
    
    for (const auto& el : this->photonDict) {
        
        Vector3d posNode = el.second;
        double d;
        d = sqrt((-posNode.x/kpc -pos.x/kpc)*(-posNode.x/kpc-pos.x/kpc)+(posNode.y/kpc-pos.y/kpc)*(posNode.y/kpc-pos.y/kpc)+(posNode.z/kpc-pos.z/kpc)*(posNode.z/kpc-pos.z/kpc));
        
        if (d<dMin) {
            dMin = d;
            iMin = el.first;
        }
    }
    
    if (iMin == -1) {
        return -1.;
    } else {
       
        if ((ePhoton < photonEnergies[iMin][0]) || (ePhoton > photonEnergies[iMin][photonEnergies[iMin].size() - 1])){
            return 0;
        } else {
            std::vector<double> rowE = this->photonEnergies[iMin];
            std::vector<double> rowD = this->photonDensity[iMin];
            return interpolate(ePhoton, rowE, rowD);
        }
    }
}

double TabularSpatialPhotonField::getMinimumPhotonEnergy(double z, const Vector3d &pos) const {
    
    double dMin = 1000.;
    int iMin;
    
    for (const auto& el : this->photonDict) {
        
        Vector3d posNode = el.second;
        double d;
        d = sqrt((-posNode.x/kpc -pos.x/kpc)*(-posNode.x/kpc-pos.x/kpc)+(posNode.y/kpc-pos.y/kpc)*(posNode.y/kpc-pos.y/kpc)+(posNode.z/kpc-pos.z/kpc)*(posNode.z/kpc-pos.z/kpc));
        
        if (d<dMin) {
            dMin = d;
            iMin = el.first;
        }
    }
    return photonEnergies[iMin][0];
}

double TabularSpatialPhotonField::getMaximumPhotonEnergy(double z, const Vector3d &pos) const {
    
    double dMin = 1000.;
    int iMin;
    
    for (const auto& el : this->photonDict) {
        
        Vector3d posNode = el.second;
        double d;
        d = sqrt((-posNode.x/kpc -pos.x/kpc)*(-posNode.x/kpc-pos.x/kpc)+(posNode.y/kpc-pos.y/kpc)*(posNode.y/kpc-pos.y/kpc)+(posNode.z/kpc-pos.z/kpc)*(posNode.z/kpc-pos.z/kpc));
        
        if (d<dMin) {
            dMin = d;
            iMin = el.first;
        }
    }
    return photonEnergies[iMin][photonEnergies[iMin].size() - 1];
}

std::vector<double> TabularSpatialPhotonField::readPhotonEnergy(std::string filePath) {
    std::ifstream infile(filePath.c_str());
    if (!infile.good())
        throw std::runtime_error("TabularPhotonField::readPhotonEnergy: could not open " + filePath);

    std::string line;
    std::vector<double> vE;
    
    while (std::getline(infile, line)) {
        if ((line.size() > 0) & (line[0] != '#') ) {
            vE.insert(vE.begin(),std::stod(line));
        }
    }
    infile.close();
    return vE;
}

std::vector<double> TabularSpatialPhotonField::readPhotonDensity(std::string filePath) {
    std::ifstream infile(filePath.c_str());
    if (!infile.good())
        throw std::runtime_error("TabularPhotonField::readPhotonDensity: could not open " + filePath);

    std::string line;
    std::vector<double> vD;
    
    while (std::getline(infile, line)) {
        if ((line.size() > 0) & (line[0] != '#') )
            vD.insert(vD.begin(),std::stod(line));
    }
    infile.close();
    return vD;
}

void TabularSpatialPhotonField::checkInputData() const {
    
    std::size_t numRowsDens = this->photonEnergies.size();
    std::size_t numRowsEn = this->photonDensity.size();
    
    if (numRowsDens == numRowsEn) {
        for (int j = 0; j < this->photonDensity.size(); ++j) { //take the proper row size!
            if (this->photonEnergies[j].size() != this->photonDensity[j].size())
                throw std::runtime_error("TabularPhotonField::checkInputData: length of photon energy input is unequal to length of photon density input");
            for (int i = 0; i < this->photonEnergies[j].size(); ++i) {
                double ePrevious = 0.;
                double e = this->photonEnergies[j][i];
                if (e <= 0.)
                    throw std::runtime_error("TabularSpatialPhotonField::checkInputData: a value in the photon energy input is not positive");
                if (e <= ePrevious)
                    throw std::runtime_error("TabularSpatialPhotonField::checkInputData: photon energy values are not strictly increasing");
                ePrevious = e;
            }
            
            for (int i = 0; i < this->photonDensity[j].size(); ++i) {
                if (this->photonDensity[j][i] < 0.)
                    throw std::runtime_error("TabularSpatialPhotonField::checkInputData: a value in the photon density input is negative");
            }
        }
    } else {
        throw std::runtime_error("TabularSpatialPhotonField::checkInputData: different number of rows in this->photonEnergies and this->photonDensity.");
    }
}

BlackbodyPhotonField::BlackbodyPhotonField(std::string fieldName, double blackbodyTemperature) {
    this->fieldName = fieldName;
    this->blackbodyTemperature = blackbodyTemperature;
    this->quantile = 0.0001; // tested to be sufficient, only used for extreme values of primary energy or temperature
}

double BlackbodyPhotonField::getPhotonDensity(double Ephoton, double z, const Vector3d &pos) const {
	return 8 * M_PI * pow_integer<3>(Ephoton / (h_planck * c_light)) / std::expm1(Ephoton / (k_boltzmann * this->blackbodyTemperature));
}

double BlackbodyPhotonField::getMinimumPhotonEnergy(double z, const Vector3d &pos) const {
	double A;
	int quantile_int = 10000 * quantile;
	switch (quantile_int)
	{
	case 1:	// 0.01 % percentil
		A = 1.093586e-5 * eV / kelvin;
		break;
	case 10:		// 0.1 % percentil
		A = 2.402189e-5 * eV / kelvin;
		break;
	case 100:		// 1 % percentil
		A = 5.417942e-5 * eV / kelvin;
		break;
	default:
		throw std::runtime_error("Quantile not understood. Please use 0.01 (1%), 0.001 (0.1%) or 0.0001 (0.01%) \n");
		break;
	}
	return A * this -> blackbodyTemperature;
}

double BlackbodyPhotonField::getMaximumPhotonEnergy(double z, const Vector3d &pos) const {
	double factor = std::max(1., blackbodyTemperature / 2.73);
	return 0.1 * factor * eV; // T dependent scaling, starting at 0.1 eV as suitable for CMB
}

void BlackbodyPhotonField::setQuantile(double q) {
	if(not ((q == 0.0001) or (q == 0.001) or (q == 0.01)))
		throw std::runtime_error("Quantile not understood. Please use 0.01 (1%), 0.001 (0.1%) or 0.0001 (0.01%) \n");
	this -> quantile = q;
}

} // namespace crpropa
