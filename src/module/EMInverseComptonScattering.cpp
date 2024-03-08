#include "crpropa/module/EMInverseComptonScattering.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/Common.h"

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

static const double mec2 = mass_electron * c_squared;

EMInverseComptonScattering::EMInverseComptonScattering(ref_ptr<PhotonField> photonField, bool havePhotons, double thinning, double limit) {
	setPhotonField(photonField);
	setHavePhotons(havePhotons);
	setLimit(limit);
	setThinning(thinning);
}

void EMInverseComptonScattering::setPhotonField(ref_ptr<PhotonField> photonField) {
	
    this->photonField = photonField;
	std::string fname = photonField->getFieldName();
	setDescription("EMInverseComptonScattering: " + fname);
    
    if (!this->photonField->hasPositionDependence()) {
        
        this->interactionRates = new InteractionRatesIsotropic("interactionRatesIsotropic", false);
        InteractionRatesIsotropic* intRatesIso = static_cast<InteractionRatesIsotropic*>(this->interactionRates.get()); //there's the dedicated function in CRPropa
        
        initRate(getDataPath("EMInverseComptonScattering/rate_" + fname + ".txt"), intRatesIso);
        initCumulativeRate(getDataPath("EMInverseComptonScattering/cdf_" + fname + ".txt"), intRatesIso);
        
    } else {
        
        this->interactionRates = new InteractionRatesPositionDependent("interactionRatesPositionDependent", true);
        InteractionRatesPositionDependent* intRatesPosDep = static_cast<InteractionRatesPositionDependent*>(this->interactionRates.get());
        
        initRatePositionDependentPhotonField(getDataPath("EMInverseComptonScattering/"+fname+"/Rate/"), intRatesPosDep);
        initCumulativeRatePositionDependentPhotonField(getDataPath("EMInverseComptonScattering/"+fname+"/CumulativeRate/"), intRatesPosDep);
        
    }
}

void EMInverseComptonScattering::setHavePhotons(bool havePhotons) {
	this->havePhotons = havePhotons;
}

void EMInverseComptonScattering::setLimit(double limit) {
	this->limit = limit;
}

void EMInverseComptonScattering::setThinning(double thinning) {
	this->thinning = thinning;
}

void EMInverseComptonScattering::initRate(std::string filename, InteractionRatesIsotropic* intRatesIso) {
	std::ifstream infile(filename.c_str());

    std::vector<double> tabEnergy;
    std::vector<double> tabRate;
    
	if (!infile.good())
		throw std::runtime_error("EMInverseComptonScattering: could not open file " + filename);

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
    
    intRatesIso->setTabulatedEnergy(tabEnergy);
    intRatesIso->setTabulatedRate(tabRate);
}

std::string EMInverseComptonScattering::splitFilename(const std::string str) {
            std::size_t found = str.find_last_of("/\\");
            std::string s = str.substr(found+1);
            return s;
    }

void EMInverseComptonScattering::initRatePositionDependentPhotonField(std::string filepath, InteractionRatesPositionDependent* intRatesPosDep) {
    
    std::vector<std::vector<double>> tabEnergy;
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
            std::runtime_error("EMInverseComptonScattering: could not open file " + filename);
        
        while (infile.good()) {
            if (infile.peek() != '#') {
                double a, b;
                infile >> a >> b;
                if (infile) {
                    vecEnergy.push_back(pow(10, a) * eV);
                    vecRate.push_back(b / Mpc);
                }
            }
            infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
        }
        
        tabEnergy.push_back(vecEnergy);
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
    
    intRatesPosDep->setTabulatedEnergy(tabEnergy);
    intRatesPosDep->setTabulatedRate(tabRate);
    intRatesPosDep->setPhotonDict(photonDict);
    
}

void EMInverseComptonScattering::initCumulativeRate(std::string filename, InteractionRatesIsotropic* intRatesIso) {
	
    std::ifstream infile(filename.c_str());

    std::vector<double> tabE;
    std::vector<double> tabs;
    std::vector<std::vector<double>> tabCDF;
    
	if (!infile.good())
		throw std::runtime_error("EMInverseComptonScattering: could not open file " + filename);
	
	// skip header
	while (infile.peek() == '#')
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');

	// read s values in first line
	double a;
	infile >> a; // skip first value
	while (infile.good() and (infile.peek() != '\n')) {
		infile >> a;
		tabs.push_back(pow(10, a) * eV * eV);
	}

	// read all following lines: E, cdf values
	while (infile.good()) {
		infile >> a;
		if (!infile)
			break;  // end of file
		tabE.push_back(pow(10, a) * eV);
		std::vector<double> cdf;
		for (int i = 0; i < tabs.size(); i++) {
			infile >> a;
			cdf.push_back(a / Mpc);
		}
		tabCDF.push_back(cdf);
	}
	infile.close();
    
    intRatesIso->setTabulatedE(tabE);
    intRatesIso->setTabulateds(tabs);
    intRatesIso->setTabulatedCDF(tabCDF);
}

void EMInverseComptonScattering::initCumulativeRatePositionDependentPhotonField(std::string filepath, InteractionRatesPositionDependent* intRatesPosDep) {
    
    std::vector<std::vector<double>> tabE;
    std::vector<std::vector<double>> tabs;
    std::vector<std::vector<std::vector<double>>> tabCDF;
    
    std::__fs::filesystem::path dir = filepath;
    int iFile = 0;
    
    for (auto const& dir_entry : std::__fs::filesystem::directory_iterator{dir}) {
        
        std::vector<double> vecE;
        std::vector<double> vecs;
        std::vector<std::vector<double>> vecCDF;
        
        // the input filename here should be a string
        //check if it is correct, i.e. a proper filename string
        std::string filename = dir_entry.path().string();
        std::ifstream infile(filename.c_str());
        
        if (!infile.good())
            throw std::runtime_error("EMInverseComptonScattering: could not open file " + filename);
        
        // skip header
        while (infile.peek() == '#')
            infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
        
        // read s values in first line
        double a;
        infile >> a; // skip first value
        while (infile.good() and (infile.peek() != '\n')) {
            infile >> a;
            vecs.push_back(pow(10, a) * eV * eV);
        }
        
        // read all following lines: E, cdf values
        while (infile.good()) {
            infile >> a;
            if (!infile)
                break;  // end of file
            vecE.push_back(pow(10, a) * eV);
            std::vector<double> cdf;
            for (int i = 0; i < tabs.size(); i++) {
                infile >> a;
                cdf.push_back(a / Mpc);
            }
            vecCDF.push_back(cdf);
        }
        iFile = iFile + 1;
        tabE.push_back(vecE);
        tabs.push_back(vecs);
        tabCDF.push_back(vecCDF);
        infile.close();
    }
    
    intRatesPosDep->setTabulatedE(tabE);
    intRatesPosDep->setTabulateds(tabs);
    intRatesPosDep->setTabulatedCDF(tabCDF);
}

// Class to calculate the energy distribution of the ICS photon and to sample from it
class ICSSecondariesEnergyDistribution {
	private:
		std::vector< std::vector<double> > data;
		std::vector<double> s_values;
		size_t Ns;
		size_t Nrer;
		double s_min;
		double s_max;
		double dls;

	public:
		// differential cross-section, see Lee '96 (arXiv:9604098), eq. 23 for x = Ee'/Ee
		double dSigmadE(double x, double beta) {
			double q = ((1 - beta) / beta) * (1 - 1./x);
			return ((1 + beta) / beta) * (x + 1./x + 2 * q + q * q);
		}

		// create the cumulative energy distribution of the up-scattered photon
		ICSSecondariesEnergyDistribution() {
			Ns = 1000;
			Nrer = 1000;
			s_min = mec2 * mec2;
			s_max = 2e23 * eV * eV;
			dls = (log(s_max) - log(s_min)) / Ns;
			data = std::vector< std::vector<double> >(1000, std::vector<double>(1000));
			std::vector<double> data_i(1000);

			// tabulate s bin borders
			s_values = std::vector<double>(1001);
			for (size_t i = 0; i < Ns + 1; ++i)
				s_values[i] = s_min * exp(i*dls);


			// for each s tabulate cumulative differential cross section
			for (size_t i = 0; i < Ns; i++) {
				double s = s_min * exp((i+0.5) * dls);
				double beta = (s - s_min) / (s + s_min);
				double x0 = (1 - beta) / (1 + beta);
				double dlx = -log(x0) / Nrer;

				// cumulative midpoint integration
				data_i[0] = dSigmadE(x0, beta) * expm1(dlx);
				for (size_t j = 1; j < Nrer; j++) {
					double x = x0 * exp((j+0.5) * dlx);
					double dx = exp((j+1) * dlx) - exp(j * dlx);
					data_i[j] = dSigmadE(x, beta) * dx;
					data_i[j] += data_i[j-1];
				}
				data[i] = data_i;
			}
		}

		// draw random energy for the up-scattered photon Ep(Ee, s)
		double sample(double Ee, double s) {
			size_t idx = std::lower_bound(s_values.begin(), s_values.end(), s) - s_values.begin();
			std::vector<double> s0 = data[idx];
			Random &random = Random::instance();
			size_t j = random.randBin(s0) + 1; // draw random bin (upper bin boundary returned)
			double beta = (s - s_min) / (s + s_min);
			double x0 = (1 - beta) / (1 + beta);
			double dlx = -log(x0) / Nrer;
			double binWidth = x0 * (exp(j * dlx) - exp((j-1) * dlx));
			double Ep = (x0 * exp((j-1) * dlx) + binWidth) * Ee;
			return std::min(Ee, Ep); // prevent Ep > Ee from numerical inaccuracies
		}
};

void EMInverseComptonScattering::getPerformInteractionTabs(const Vector3d &position, std::vector<double> &tabE, std::vector<double> &tabs, std::vector<std::vector<double>> &tabCDF) const {
    if (!this->photonField->hasPositionDependence()){
        
        InteractionRatesIsotropic* intRateIso = static_cast<InteractionRatesIsotropic*>(this->interactionRates.get());
        
        tabE = intRateIso->getTabulatedE();
        tabs = intRateIso->getTabulateds();
        tabCDF = intRateIso->getTabulatedCDF();
        
    } else {
        
        InteractionRatesPositionDependent* intRatePosDep = static_cast<InteractionRatesPositionDependent*>(this->interactionRates.get());
        
        std::vector<std::vector<double>> E = intRatePosDep->getTabulatedE();
        std::vector<std::vector<double>> s = intRatePosDep->getTabulateds();
        std::vector<std::vector<std::vector<double>>> CDF = intRatePosDep->getTabulatedCDF();
        std::unordered_map<int,Vector3d> photonDict = intRatePosDep->getPhotonDict();
        
        double dMin = 1000. * kpc;
        int iMin = -1;
        
        for (const auto& el : photonDict) {
            
            Vector3d posNode = el.second;
            double d;
            d = sqrt((- posNode.x / kpc - position.x / kpc) * (- posNode.x / kpc - position.x / kpc) + (posNode.y / kpc - position.y / kpc) * (posNode.y / kpc - position.y / kpc) + (posNode.z / kpc - position.z / kpc) * (posNode.z / kpc - position.z / kpc));
            
            if (d<dMin) {
                dMin = d;
                iMin = el.first;
            }
        }
        
        tabE = E[iMin];
        tabs = s[iMin];
        tabCDF = CDF[iMin];
    }
}

void EMInverseComptonScattering::getProcessTabs(const Vector3d &position, std::vector<double> &tabEnergy, std::vector<double> &tabRate) const {
    if (!this->photonField->hasPositionDependence()) {
        
        InteractionRatesIsotropic* intRateIso = static_cast<InteractionRatesIsotropic*>(this->interactionRates.get());
        
        tabEnergy = intRateIso->getTabulatedEnergy();
        tabRate = intRateIso->getTabulatedRate();
        
    } else {
        
        InteractionRatesPositionDependent* intRatePosDep = static_cast<InteractionRatesPositionDependent*>(this->interactionRates.get());
        
        std::vector<std::vector<double>> Energy = intRatePosDep->getTabulatedEnergy();
        std::vector<std::vector<double>> Rate = intRatePosDep->getTabulatedRate();
        std::unordered_map<int,Vector3d> photonDict = intRatePosDep->getPhotonDict();
        
        double dMin = 1000. * kpc;
        int iMin = -1;
        
        for (const auto& el : photonDict) {
            
            Vector3d posNode = el.second;
            double d;
            d = sqrt((- posNode.x / kpc - position.x / kpc) * (- posNode.x / kpc - position.x / kpc) + (posNode.y / kpc - position.y / kpc) * (posNode.y / kpc - position.y / kpc) + (posNode.z / kpc - position.z / kpc) * (posNode.z / kpc - position.z / kpc));
            
            if (d<dMin) {
                dMin = d;
                iMin = el.first;
            }
        }
        
        tabEnergy = Energy[iMin];
        tabRate = Rate[iMin];
    }
}

void EMInverseComptonScattering::performInteraction(Candidate *candidate) const {
	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

    Vector3d position = candidate->current.getPosition();
    
    std::vector<double> tabE;
    std::vector<double> tabs;
    std::vector<std::vector<double>> tabCDF;
    
    getPerformInteractionTabs(position, tabE, tabs, tabCDF);
    
	if (E < tabE.front() or E > tabE.back())
		return;

	// sample the value of s
	Random &random = Random::instance();
	size_t i = closestIndex(E, tabE);
	size_t j = random.randBin(tabCDF[i]);
	double s_kin = pow(10, log10(tabs[j]) + (random.rand() - 0.5) * 0.1);
	double s = s_kin + mec2 * mec2;

	// sample electron energy after scattering
	static ICSSecondariesEnergyDistribution distribution;
	double Enew = distribution.sample(E, s);

	// add up-scattered photon
	double Esecondary = E - Enew;
	double f = Enew / E;
	if (havePhotons) {
		if (random.rand() < pow(1 - f, thinning)) {
			double w = 1. / pow(1 - f, thinning);
			Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
			candidate->addSecondary(22, Esecondary / (1 + z), pos, w, interactionTag);
		}
	}

	// update the primary particle energy; do this after adding the secondary to correctly set the secondary's parent
	candidate->current.setEnergy(Enew / (1 + z));
}

void EMInverseComptonScattering::process(Candidate *candidate) const {
	// check if electron / positron
	int id = candidate->current.getId();
	if (abs(id) != 11)
		return;

	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);
    Vector3d position = candidate->current.getPosition();

    std::vector<double> tabEnergy;
    std::vector<double> tabRate;
    
    getProcessTabs(position, tabEnergy, tabRate);

	if (E < tabEnergy.front() or (E > tabEnergy.back()))
		return;

	// interaction rate
	double rate = interpolate(E, tabEnergy, tabRate);
	rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);

	// run this loop at least once to limit the step size
	double step = candidate->getCurrentStep();
	Random &random = Random::instance();
	do {
		double randDistance = -log(random.rand()) / rate;

		// check for interaction; if it doesn't ocurr, limit next step
		if (step < randDistance) {
			candidate->limitNextStep(limit / rate);
			return;
		}
		performInteraction(candidate);

		// repeat with remaining step
		step -= randDistance;
	} while (step > 0);
}

void EMInverseComptonScattering::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string EMInverseComptonScattering::getInteractionTag() const {
	return interactionTag;
}

} // namespace crpropa
