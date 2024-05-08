#include "crpropa/module/EMDoublePairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/InteractionRates.h"

#include <fstream>
#include <limits>
#include <stdexcept>
#include <filesystem>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace crpropa {

EMDoublePairProduction::EMDoublePairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons, double thinning, double limit) {
	setPhotonField(photonField);
	setHaveElectrons(haveElectrons);
	setLimit(limit);
	setThinning(thinning);
}

void EMDoublePairProduction::setPhotonField(ref_ptr<PhotonField> photonField) {
	
    this->photonField = photonField;
	std::string fname = photonField->getFieldName();
	setDescription("EMDoublePairProduction: " + fname);
    if (!this->photonField->hasPositionDependence()) {
        
        this->interactionRates = new InteractionRatesHomogeneous("interactionRatesHomogeneous", false);
        InteractionRatesHomogeneous* intRatesHom = static_cast<InteractionRatesHomogeneous*>(this->interactionRates.get());
        initRate(getDataPath("EMDoublePairProduction/rate_" + fname + ".txt"), intRatesHom);
        
    } else {
        
        this->interactionRates = new InteractionRatesPositionDependent("interactionRatesPositionDependent", true);
        InteractionRatesPositionDependent* intRatesPosDep = static_cast<InteractionRatesPositionDependent*>(this->interactionRates.get());
        initRatePositionDependentPhotonField(getDataPath("EMDoublePairProduction/"+fname+"/Rate/"), intRatesPosDep);
        
    }
}

void EMDoublePairProduction::setHaveElectrons(bool haveElectrons) {
	this->haveElectrons = haveElectrons;
}

void EMDoublePairProduction::setLimit(double limit) {
	this->limit = limit;
}

void EMDoublePairProduction::setThinning(double thinning) {
	this->thinning = thinning;
}

void EMDoublePairProduction::initRate(std::string filename, InteractionRatesHomogeneous* intRatesHom) {
	std::ifstream infile(filename.c_str());

    std::vector<double> tabEnergy;
    std::vector<double> tabRate;
    
	if (!infile.good())
		throw std::runtime_error("EMDoublePairProduction: could not open file " + filename);

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabEnergy.push_back(pow(10, a) * eV);
				tabRate.push_back(b / Mpc);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
    
    intRatesHom->setTabulatedEnergy(tabEnergy);
    intRatesHom->setTabulatedRate(tabRate);
    
}

std::string EMDoublePairProduction::splitFilename(const std::string str) {
            std::size_t found = str.find_last_of("/\\");
            std::string s = str.substr(found+1);
            return s;
}

void EMDoublePairProduction::initRatePositionDependentPhotonField(std::string filepath, InteractionRatesPositionDependent* intRatesPosDep) {

    std::vector<std::vector<double>> tabRate;
    
    std::__fs::filesystem::path dir = filepath;
    std::unordered_map<int, Vector3d> photonDict;
    int iFile = 0;
    
    for (auto const& dir_entry : std::__fs::filesystem::directory_iterator{dir}) {

        // the input filename here should be a string
        //check if it is correct, i.e. a proper filename string
        std::string filename = dir_entry.path().string();
        std::ifstream infile(filename.c_str());
        
        std::vector<double> vecEnergy;
        std::vector<double> vecRate;
        
        if (!infile.good())
            throw
            std::runtime_error("EMDoublePairProduction: could not open file " + filename);
        
        while (infile.good()) {
            if (infile.peek() != '#') {
                double a, b;
                infile >> a >> b;
                if (infile) {
                    if (iFile == 0) {
                        vecEnergy.push_back(pow(10, a) * eV);
                        intRatesPosDep->setTabulatedEnergy(vecEnergy);
                    }
                    vecRate.push_back(b / Mpc);
                }
            }
            infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
        }
        
        tabRate.push_back(vecRate);
        
        double x, y, z;
        std::string str;
        std::stringstream ss;
        
        std::string filename_split = splitFilename(dir_entry.path().string());
        ss << filename_split;
        
        int iLine = 0;
        
        while (getline(ss, str, '_')) {
            if (iLine == 3) {
                x = stod(str) * kpc;
            }
            if (iLine == 4) {
                y = stod(str) * kpc;
            }
            if (iLine == 5) {
                z = stod(str) * kpc;
            }
            iLine = iLine + 1;
        }
        
        Vector3d vPos(x, y, z);
        photonDict[iFile] = vPos;
        
        iFile = iFile + 1;
        infile.close();
    }

    intRatesPosDep->setTabulatedRate(tabRate);
    intRatesPosDep->setPhotonDict(photonDict);
    
}

void EMDoublePairProduction::getProcessTabs(const Vector3d &position, std::vector<double> &tabEnergy, std::vector<double> &tabRate) const {
    if (!this->photonField->hasPositionDependence()) {
        
        InteractionRatesHomogeneous* intRateHom = static_cast<InteractionRatesHomogeneous*>(this->interactionRates.get());
        
        tabEnergy = intRateHom->getTabulatedEnergy();
        tabRate = intRateHom->getTabulatedRate();
        
    } else {
        
        InteractionRatesPositionDependent* intRatePosDep = static_cast<InteractionRatesPositionDependent*>(this->interactionRates.get());
        
        std::vector<double> Energy = intRatePosDep->getTabulatedEnergy();
        std::vector<std::vector<double>> Rate = intRatePosDep->getTabulatedRate();
        std::unordered_map<int,Vector3d> photonDict = intRatePosDep->getPhotonDict();
        
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
        
        tabEnergy = Energy;
        tabRate = Rate[iMin];
    }
}

void EMDoublePairProduction::performInteraction(Candidate *candidate) const {
	// the photon is lost after the interaction
	candidate->setActive(false);

	if (not haveElectrons)
		return;

	// Use assumption of Lee 96 arXiv:9604098
	// Energy is equally shared between one e+e- pair, but take mass of second e+e- pair into account.
	// This approximation has been shown to be valid within -1.5%.
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);
	double Ee = (E - 2 * mass_electron * c_squared) / 2;

	Random &random = Random::instance();
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

	double f = Ee / E;

	if (haveElectrons) {
		if (random.rand() < pow(1 - f, thinning)) {
			double w = 1. / pow(1 - f, thinning);
			candidate->addSecondary( 11, Ee / (1 + z), pos, w, interactionTag);
		} 
		if (random.rand() < pow(f, thinning)) {
			double w = 1. / pow(f, thinning);
			candidate->addSecondary(-11, Ee / (1 + z), pos, w, interactionTag);
		}
	}
}

void EMDoublePairProduction::process(Candidate *candidate) const {
	// check if photon
	if (candidate->current.getId() != 22)
		return;

	// scale the electron energy instead of background photons
	double z = candidate->getRedshift();
	double E = (1 + z) * candidate->current.getEnergy();
    Vector3d position = candidate->current.getPosition();

    std::vector<double> tabEnergy;
    std::vector<double> tabRate;
    
    getProcessTabs(position, tabEnergy, tabRate);
    
	// check if in tabulated energy range
	if (E < tabEnergy.front() or (E > tabEnergy.back()))
		return;

	// interaction rate
	double rate = interpolate(E, tabEnergy, tabRate);
	rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);

	// check for interaction
	Random &random = Random::instance();
	double randDistance = -log(random.rand()) / rate;
	double step = candidate->getCurrentStep();
	if (step < randDistance) {
		candidate->limitNextStep(limit / rate);
		return;
	} else { // after performing interaction photon ceases to exist (hence return)
		performInteraction(candidate);
		return;
	}

}

void EMDoublePairProduction::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string EMDoublePairProduction::getInteractionTag() const {
	return interactionTag;
}


} // namespace crpropa
