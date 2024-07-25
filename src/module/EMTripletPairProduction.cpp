#include "crpropa/module/EMTripletPairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

static const double mec2 = mass_electron * c_squared;

EMTripletPairProduction::EMTripletPairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons, double thinning, double limit) {
	setPhotonField(photonField);
	setHaveElectrons(haveElectrons);
	setLimit(limit);
	setThinning(thinning);
}

void EMTripletPairProduction::setPhotonField(ref_ptr<PhotonField> photonField) {
	this->photonField = photonField;
	std::string fname = photonField->getFieldName();
	setDescription("EMTripletPairProduction: " + fname);
    if (!this->photonField->hasPositionDependence()){
        
        this->interactionRates = new InteractionRatesHomogeneous("interactionRatesHomogeneous", false);
        InteractionRatesHomogeneous* intRatesHom = static_cast<InteractionRatesHomogeneous*>(this->interactionRates.get());
        
        initRate(getDataPath("EMTripletPairProduction/rate_" + fname + ".txt"), intRatesHom);
        initCumulativeRate(getDataPath("EMTripletPairProduction/cdf_" + fname + ".txt"), intRatesHom);
        
    } else {
        
        this->interactionRates = new InteractionRatesPositionDependent("interactionRatesPositionDependent", true);
        InteractionRatesPositionDependent* intRatesPosDep = static_cast<InteractionRatesPositionDependent*>(this->interactionRates.get());
        
        initRatePositionDependentPhotonField(getDataPath("EMTripletPairProduction/"+fname+"/Rate/"), intRatesPosDep);
        initCumulativeRatePositionDependentPhotonField(getDataPath("EMTripletPairProduction/"+fname+"/CumulativeRate/"), intRatesPosDep);
        
    }
}

void EMTripletPairProduction::setHaveElectrons(bool haveElectrons) {
	this->haveElectrons = haveElectrons;
}

void EMTripletPairProduction::setLimit(double limit) {
	this->limit = limit;
}

void EMTripletPairProduction::setThinning(double thinning) {
	this->thinning = thinning;
}

void EMTripletPairProduction::initRate(std::string filename, InteractionRatesHomogeneous* intRatesHom) {
	std::ifstream infile(filename.c_str());

    std::vector<double> tabEnergy;
    std::vector<double> tabRate;
    
	if (!infile.good())
		throw std::runtime_error("EMTripletPairProduction: could not open file " + filename);

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

std::string EMTripletPairProduction::splitFilename(const std::string str) {
            std::size_t found = str.find_last_of("/\\");
            std::string s = str.substr(found+1);
            return s;
}

void EMTripletPairProduction::initRatePositionDependentPhotonField(std::string filepath, InteractionRatesPositionDependent* intRatesPosDep) {
    
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
            std::runtime_error("EMTripletPairProduction: could not open file " + filename);
        
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
     
void EMTripletPairProduction::initCumulativeRate(std::string filename, InteractionRatesHomogeneous* intRatesHom) {
	std::ifstream infile(filename.c_str());

    std::vector<double> tabE;
    std::vector<double> tabs;
    std::vector<std::vector<double>> tabCDF;
    
	if (!infile.good())
		throw std::runtime_error(
				"EMTripletPairProduction: could not open file " + filename);

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

void EMTripletPairProduction::initCumulativeRatePositionDependentPhotonField(std::string filepath, InteractionRatesPositionDependent* intRatesPosDep) {
    
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
            throw std::runtime_error("EMTripletPairProduction: could not open file " + filename);
        
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
 
    intRatesPosDep->setTabulateds(tabs);
    intRatesPosDep->setTabulatedCDF(tabCDF);
}

void EMTripletPairProduction::performInteraction(Candidate *candidate) const {
	int id = candidate->current.getId();
	if  (abs(id) != 11)
		return;

	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);
    Vector3d position = candidate->current.getPosition();
    
    std::vector<double> tabE;
    std::vector<double> tabs;
    std::vector<std::vector<double>> tabCDF;
    
    this->interactionRates->loadPerformInteractionTabs(position, tabE, tabs, tabCDF);
	
    if (E < tabE.front() or E > tabE.back())
		return;

	// sample the value of eps
	Random &random = Random::instance();
	size_t i = closestIndex(E, tabE);
	size_t j = random.randBin(tabCDF[i]);
	double s_kin = pow(10, log10(tabs[j]) + (random.rand() - 0.5) * 0.1);
	double eps = s_kin / 4. / E; // random background photon energy

	// Use approximation from A. Mastichiadis et al., Astroph. Journ. 300:178-189 (1986), eq. 30.
	// This approx is valid only for alpha >=100 where alpha = p0*eps*costheta - E0*eps
	// For our purposes, me << E0 --> p0~E0 --> alpha = E0*eps*(costheta - 1) >= 100
	double Epp = 5.7e-1 * pow(eps / mec2, -0.56) * pow(E / mec2, 0.44) * mec2;

	double f = Epp / E;

	if (haveElectrons) {
		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
		if (random.rand() < pow(1 - f, thinning)) {
			double w = 1. / pow(1 - f, thinning);
			candidate->addSecondary(11, Epp / (1 + z), pos, w, interactionTag);
		}
		if (random.rand() < pow(f, thinning)) {
			double w = 1. / pow(f, thinning);
			candidate->addSecondary(-11, Epp / (1 + z), pos, w, interactionTag);
		}
	}
	// Update the primary particle energy.
	// This is done after adding the secondaries to correctly set the secondaries parent
	candidate->current.setEnergy((E - 2 * Epp) / (1. + z));
}

void EMTripletPairProduction::process(Candidate *candidate) const {
	// check if electron / positron
	int id = candidate->current.getId();
	if (abs(id) != 11)
		return;

	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = (1 + z) * candidate->current.getEnergy();
    Vector3d position = candidate->current.getPosition();
    
    double scaling = pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);
    double rate = this->interactionRates->getProcessRate(E, position);
    
    rate *= scaling;

	// run this loop at least once to limit the step size
	double step = candidate->getCurrentStep();
	Random &random = Random::instance();
	do {
		double randDistance = -log(random.rand()) / rate;
		// check for interaction; if it doesn't occur, limit next step
		if (step < randDistance) { 
			candidate->limitNextStep(limit / rate);
			return;
		}
		performInteraction(candidate);
		step -= randDistance; 
	} while (step > 0.);
}

void EMTripletPairProduction::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string EMTripletPairProduction::getInteractionTag() const {
	return interactionTag;
}

} // namespace crpropa
