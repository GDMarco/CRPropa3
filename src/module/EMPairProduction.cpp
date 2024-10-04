#include "crpropa/module/EMPairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/InteractionRates.h"
#include "crpropa/Geometry.h"

#include <fstream>
#include <locale>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <filesystem>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>


namespace crpropa {

static const double mec2 = mass_electron * c_squared;

EMPairProduction::EMPairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons, double thinning, double limit, Surface* surface) {
    setSurface(surface);
    setPhotonField(photonField);
	setThinning(thinning);
	setLimit(limit);
	setHaveElectrons(haveElectrons);
}

void EMPairProduction::setPhotonField(ref_ptr<PhotonField> photonField) {
    
	this->photonField = photonField;
	std::string fname = photonField->getFieldName();
	setDescription("EMPairProduction: " + fname);
    
    if (!this->photonField->hasPositionDependence()){
        
        this->interactionRates = new InteractionRatesHomogeneous("interactionRatesHomogeneous", false);
        InteractionRatesHomogeneous* intRatesHom = static_cast<InteractionRatesHomogeneous*>(this->interactionRates.get());
        
        initRate(getDataPath("EMPairProduction/rate_" + fname + ".txt"), intRatesHom);
        initCumulativeRate(getDataPath("EMPairProduction/cdf_" + fname + ".txt"), intRatesHom);
        
    } else {
        
        this->interactionRates = new InteractionRatesPositionDependent("interactionRatesPositionDependent", true);
        InteractionRatesPositionDependent* intRatesPosDep = static_cast<InteractionRatesPositionDependent*>(this->interactionRates.get());
        
        initRatePositionDependentPhotonField(getDataPath("EMPairProduction/"+fname+"/Rate/"), intRatesPosDep);
        initCumulativeRatePositionDependentPhotonField(getDataPath("EMPairProduction/"+fname+"/CumulativeRate/"), intRatesPosDep);
        
    }
}

void EMPairProduction::setHaveElectrons(bool haveElectrons) {
	this->haveElectrons = haveElectrons;
}

void EMPairProduction::setLimit(double limit) {
	this->limit = limit;
}

void EMPairProduction::setThinning(double thinning) {
	this->thinning = thinning;
}

void EMPairProduction::setSurface(Surface* surface) {
    this->surface = surface;
}

bool EMPairProduction::hasSurface() const {
    return this->surface != nullptr;
}

void EMPairProduction::initRate(std::string filename, InteractionRatesHomogeneous* intRatesHom) {
	std::ifstream infile(filename.c_str());

  std::vector<double> tabEnergy;
  std::vector<double> tabRate;
    
	if (!infile.good())
		throw std::runtime_error("EMPairProduction: could not open file " + filename);

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

std::string EMPairProduction::splitFilename(const std::string str) {
            std::size_t found = str.find_last_of("/\\");
            std::string s = str.substr(found+1);
            return s;
}
                                                       
void EMPairProduction::initRatePositionDependentPhotonField(std::string filepath, InteractionRatesPositionDependent* intRatesPosDep) {
    
    std::vector<std::vector<double>> tabRate;
    
    std::__fs::filesystem::path dir = filepath;
    std::unordered_map<int, Vector3d> photonDict;
    int iFile = 0;
    int nRSups = 0;
    
    for (auto const& dir_entry : std::__fs::filesystem::directory_iterator{dir}) {

        std::string filename = dir_entry.path().string();
        std::ifstream infile(filename.c_str());
        
        if (!infile.good())
            throw
            std::runtime_error("EMPairProduction: could not open file " + filename);
        
        double x, y, z;
        std::string str;
        std::stringstream ss;
        
        std::string filename_split = splitFilename(dir_entry.path().string());
        ss << filename_split;
        
        int iLine = 0;
        
        std::locale::global(std::locale("C"));

        while (getline(ss, str, '_')) {
            if (iLine == 3) {
                x = -std::stod(str) * kpc;
            }
            if (iLine == 4) {
                y = std::stod(str) * kpc;
            }
            if (iLine == 5) {
                z = std::stod(str) * kpc;
            }
            iLine = iLine + 1;
        }
        
        Vector3d vPos(x, y, z);
        
        if (hasSurface() and !surface->isInside(vPos))
            continue;
    
        photonDict[iFile] = vPos;
        
        std::vector<double> vecEnergy;
        std::vector<double> vecRate;
        
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
        
        std::vector<double> vecEnergy1 = intRatesPosDep->getTabulatedEnergy(); //std::defaultfloat
        if (vecRate.size() != vecEnergy1.size()) {
            nRSups = nRSups + 1;
            std::cout << std::fixed << std::setprecision(7);
            std::cout << "suspVec: " << vPos.x / kpc << " " << vPos.y / kpc << " " << vPos.z / kpc << std::endl;
        }
        
        tabRate.push_back(vecRate);
        
        iFile = iFile + 1;
        infile.close();
    }
    std::cout << "# sousp rate: " << nRSups << std::endl;
    intRatesPosDep->setTabulatedRate(tabRate);
    intRatesPosDep->setPhotonDict(photonDict);
}

void EMPairProduction::initCumulativeRate(std::string filename, InteractionRatesHomogeneous* intRatesHom) {
	
  std::ifstream infile(filename.c_str());

  std::vector<double> tabE;
  std::vector<double> tabs;
  std::vector<std::vector<double>> tabCDF;
    
	if (!infile.good())
		throw std::runtime_error("EMPairProduction: could not open file " + filename);
    
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
    
  intRatesHom->setTabulatedE(tabE);
  intRatesHom->setTabulateds(tabs);
  intRatesHom->setTabulatedCDF(tabCDF);
}

void EMPairProduction::initCumulativeRatePositionDependentPhotonField(std::string filepath, InteractionRatesPositionDependent* intRatesPosDep) {
    
    std::vector<double> tabE;
    std::vector<std::vector<double>> tabs;
    std::vector<std::vector<std::vector<double>>> tabCDF;
    
    std::__fs::filesystem::path dir = filepath;
    int iFile = 0;
    int nCDFSusp = 0;

    for (auto const& dir_entry : std::__fs::filesystem::directory_iterator{dir}) {
        
        std::vector<double> vecE;
        std::vector<double> vecs;
        std::vector<std::vector<double>> vecCDF;
        
        // the input filename here should be a string
        //check if it is correct, i.e. a proper filename string
        std::string filename = dir_entry.path().string();
        std::ifstream infile(filename.c_str());
        
        if (!infile.good())
            throw std::runtime_error("EMPairProduction: could not open file " + filename);
        
        double x, y, z;
        std::string str;
        std::stringstream ss;
        
        std::string filename_split = splitFilename(dir_entry.path().string());
        ss << filename_split;
        
        int iLine = 0;
        
        std::locale::global(std::locale("C"));

        while (getline(ss, str, '_')) {
            if (iLine == 3) {
                x = -std::stod(str) * kpc;
            }
            if (iLine == 4) {
                y = std::stod(str) * kpc;
            }
            if (iLine == 5) {
                z = std::stod(str) * kpc;
            }
            iLine = iLine + 1;
        }
        
        Vector3d vPos(x, y, z);
        
        if (hasSurface() and !surface->isInside(vPos))
            continue;
        
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
            if (iFile == 0) {
                vecE.push_back(pow(10, a) * eV);
                intRatesPosDep->setTabulatedE(vecE);
            }
            std::vector<double> cdf;
            for (int i = 0; i < tabs.size(); i++) {
                infile >> a;
                cdf.push_back(a / Mpc);
            }
            vecCDF.push_back(cdf);
        }

        iFile = iFile + 1;
        
        tabs.push_back(vecs);
        tabCDF.push_back(vecCDF);
        infile.close();
    }
     std::cout << "# suspicius CDF: " << nCDFSusp << std::endl;
    intRatesPosDep->setTabulateds(tabs);
    intRatesPosDep->setTabulatedCDF(tabCDF);
}
                                                       
// Hold an data array to interpolate the energy distribution on
class PPSecondariesEnergyDistribution {
	private:
		std::vector<double> tab_s;
		std::vector< std::vector<double> > data;
		size_t N;

	public:
		// differential cross section for pair production for x = Epositron/Egamma, compare Lee 96 arXiv:9604098
		double dSigmadE_PPx(double x, double beta) {
			double A = (x / (1. - x) + (1. - x) / x );
			double B =  (1. / x + 1. / (1. - x) );
			double y = (1 - beta * beta);
			return A + y * B - y * y / 4 * B * B;
		}

		PPSecondariesEnergyDistribution() {
			N = 1000;
			size_t Ns = 1000;
			double s_min = 4 * mec2 * mec2;
			double s_max = 1e23 * eV * eV;
			double dls = log(s_max / s_min) / Ns;
			data = std::vector< std::vector<double> >(Ns, std::vector<double>(N));
			tab_s = std::vector<double>(Ns + 1);

			for (size_t i = 0; i < Ns + 1; ++i)
				tab_s[i] = s_min * exp(i*dls); // tabulate s bin borders

			for (size_t i = 0; i < Ns; i++) {
				double s = s_min * exp(i*dls + 0.5*dls);
				double beta = sqrt(1 - s_min/s);
				double x0 = (1 - beta) / 2;
				double dx = log((1 + beta) / (1 - beta)) / N;

				// cumulative midpoint integration
				std::vector<double> data_i(1000);
				data_i[0] = dSigmadE_PPx(x0, beta) * expm1(dx);
				for (size_t j = 1; j < N; j++) {
					double x = x0 * exp(j*dx + 0.5*dx);
					double binWidth = exp((j+1)*dx)-exp(j*dx);
					data_i[j] = dSigmadE_PPx(x, beta) * binWidth + data_i[j-1];
				}
				data[i] = data_i;
			}
		}

		// sample positron energy from cdf(E, s_kin)
		double sample(double E0, double s) {
			// get distribution for given s
			size_t idx = std::lower_bound(tab_s.begin(), tab_s.end(), s) - tab_s.begin();
			if (idx > data.size())
				return NAN;
				
			std::vector<double> s0 = data[idx];

			// draw random bin
			Random &random = Random::instance();
			size_t j = random.randBin(s0) + 1;

			double s_min = 4. * mec2 * mec2;
			double beta = sqrtl(1. - s_min / s);
			double x0 = (1. - beta) / 2.;
			double dx = log((1 + beta) / (1 - beta)) / N;
			double binWidth = x0 * (exp(j*dx) - exp((j-1)*dx));
			if (random.rand() < 0.5)
				return E0 * (x0 * exp((j-1) * dx) + binWidth);
			else
				return E0 * (1 - (x0 * exp((j-1) * dx) + binWidth));
		}
};

void EMPairProduction::performInteraction(Candidate *candidate) const {
    
    // scale particle energy instead of background photon energy
    double z = candidate->getRedshift();
    double E = candidate->current.getEnergy() * (1 + z);
    Vector3d position = candidate->current.getPosition();
    
    // cosmic ray photon is lost after interacting
    candidate->setActive(false);

    // check if secondary electron pair needs to be produced
    if (not haveElectrons)
        return;
    
    std::vector<double> tabE;
    std::vector<double> tabs;
    std::vector<std::vector<double>> tabCDF;
    
    this->interactionRates->loadPerformInteractionTabs(position, tabE, tabs, tabCDF);

    // check if in tabulated energy range
    if (E < tabE.front() or (E > tabE.back()))
        return;

    // sample the value of s
    Random &random = Random::instance();
    size_t i = closestIndex(E, tabE);  // find closest tabulation point
    size_t j = random.randBin(tabCDF[i]);
    double lo = std::max(4 * mec2 * mec2, tabs[j-1]);  // first s-tabulation point below min(s_kin) = (2 me c^2)^2; ensure physical value
    double hi = tabs[j];
    double s = lo + random.rand() * (hi - lo);
    
    // sample electron / positron energy
    static PPSecondariesEnergyDistribution interpolation;
    double Ee = interpolation.sample(E, s);
    double Ep = E - Ee;
    double f = Ep / E;

    // for some backgrounds Ee=nan due to precision limitations.
    if (not std::isfinite(Ee) || not std::isfinite(Ep))
        return;

	// photon is lost after interacting
	candidate->setActive(false);

	// sample random position along current step
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
	// apply sampling
	if (random.rand() < pow(f, thinning)) {
		double w = 1. / pow(f, thinning);
		candidate->addSecondary(11, Ep / (1 + z), pos, w, interactionTag);
	}
	if (random.rand() < pow(1 - f, thinning)) {
		double w = 1. / pow(1 - f, thinning);
		candidate->addSecondary(-11, Ee / (1 + z), pos, w, interactionTag);	
	}
}

void EMPairProduction::process(Candidate *candidate) const {
    
    // check if photon
    if (candidate->current.getId() != 22)
        return;

    // scale particle energy instead of background photon energy
    double z = candidate->getRedshift();
    double E = candidate->current.getEnergy() * (1 + z);
    Vector3d position = candidate->current.getPosition();
    
    // New optimization
    double rate = this->interactionRates->getProcessRate(E, position);
     
    rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);
    
    // run this loop at least once to limit the step size
    double step = candidate->getCurrentStep();
    Random &random = Random::instance();
    do {
        double randDistance = -log(random.rand()) / rate;
        // check for interaction; if it doesn't ocurr, limit next step
        if (step < randDistance) {
            candidate->limitNextStep(limit / rate);
        } else {
            performInteraction(candidate);
            return;
        }
        step -= randDistance;
    } while (step > 0.);
}

void EMPairProduction::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string EMPairProduction::getInteractionTag() const {
	return interactionTag;
}

} // namespace crpropa
