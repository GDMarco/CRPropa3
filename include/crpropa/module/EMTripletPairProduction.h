#ifndef CRPROPA_EMTRIPLETPAIRPRODUCTION_H
#define CRPROPA_EMTRIPLETPAIRPRODUCTION_H

#include <fstream>
#include <cmath>

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/InteractionRates.h"

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class EMTripletPairProduction
 @brief Electron triplet pair production of electrons with background photons.

 This module simulates electron triplet pair production of electrons with background photons for several photon fields.
 The secondary electrons from this interaction are optionally created (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 Thinning is available. A thinning of 0 means that all particles are tracked. 
 For the maximum thinning of 1, only a few representative particles are added to the list of secondaries.
 Note that for thinning>0 the output must contain the column "weights", which should be included in the post-processing.
*/
class EMTripletPairProduction: public Module {
private:
	ref_ptr<PhotonField> photonField;
	bool haveElectrons;
	double limit;
	double thinning;
	std::string interactionTag = "EMTP";
    ref_ptr<InteractionRates> interactionRates;

public:
	/** Constructor
	 @param photonField		target photon field
	 @param haveElectrons	if true, add secondary electrons as candidates
	 @param thinning		weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
	 @param limit			step size limit as fraction of mean free path
	 */
	EMTripletPairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons = false, double thinning = 0, double limit = 0.1);

	// set the target photon field
	void setPhotonField(ref_ptr<PhotonField> photonField);

	// decide if secondary electrons are added to the simulation	
	void setHaveElectrons(bool haveElectrons);

	/** limit the step to a fraction of the mean free path
	 @param limit	fraction of the mean free path, should be between 0 and 1
	*/
	void setLimit(double limit);

	/** Apply thinning with a given thinning factor
	 * @param thinning factor of thinning (0: no thinning, 1: maximum thinning)
	 */
	void setThinning(double thinning);

	/** set a custom interaction tag to trace back this interaction
	 * @param tag string that will be added to the candidate and output
	 */
	void setInteractionTag(std::string tag);
	std::string getInteractionTag() const;
	
	void initRate(std::string filename, InteractionRatesHomogeneous* intRatesHom);
	void initCumulativeRate(std::string filename, InteractionRatesHomogeneous* intRatesHom);
    
    void initRatePositionDependentPhotonField(std::string filepath, InteractionRatesPositionDependent* intRatesPosDep);
    void initCumulativeRatePositionDependentPhotonField(std::string filepath, InteractionRatesPositionDependent* intRatesPosDep);

    /**
    void getPerformInteractionTabs(const Vector3d &position, std::vector<double> &tabE, std::vector<double> &tabs, std::vector<std::vector<double>> &tabCDF) const;
    void getProcessTabs(const Vector3d &position, std::vector<double> &tabEnergy, std::vector<double> &tabRate) const;
    */
    
	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate) const;

protected:
    std::string splitFilename(const std::string str);

};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_EMTRIPLETPAIRPRODUCTION_H
