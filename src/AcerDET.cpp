#include "AcerDET.h"
#include <cstdio>

using namespace AcerDet;

/*
 * Construct AcerDET object from given configuration
 */
AcerDET::AcerDET(
	const conf::Configuration& config,
	core::IParticleDataProviderFactory *partFactory,
	core::IHistogramManager *histoManager ) :
		histos				( histoManager ),
		histos_initialized	( false )
{
	partProvider = partFactory->create();
	analyse_Cell =
		new analyse::Cell( config, histoManager, partProvider );
	analyse_BJet =
		new analyse::BJet( config, histoManager );
	analyse_Calibration =
		new analyse::Calibration( config, histoManager );
	analyse_CJet =
		new analyse::CJet( config, histoManager );
	analyse_Cluster =
		new analyse::Cluster( config, histoManager, partProvider );
	analyse_Electron =
		new analyse::Electron( config, histoManager );
	analyse_Jet =
		new analyse::Jet( config, histoManager, partProvider );
	analyse_Mis	=
		new analyse::Mis( config, histoManager );
	analyse_Muon =
		new analyse::Muon( config, histoManager );
	analyse_Photon =
		new analyse::Photon( config, histoManager );
	analyse_Tau =
		new analyse::Tau( config, histoManager );
}

/*
 * Destructor - not used
 */
AcerDET::~AcerDET() {
	histos = NULL;
	
	delete analyse_Cell;
	delete analyse_BJet;
	delete analyse_Calibration;
	delete analyse_CJet;
	delete analyse_Cluster;
	delete analyse_Electron;
	delete analyse_Jet;
	delete analyse_Mis;
	delete analyse_Muon;
	delete analyse_Photon;
	delete analyse_Tau;

}

/*
 * Analyse single InputRecord from event
 */
void AcerDET::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord, Real64_t weigth ) {

	if (!histos_initialized) {
		histos->init();
		histos_initialized = true;
	}
	
	analyse_Cell		->analyseRecord( irecord, orecord, weigth );
	analyse_Cluster		->analyseRecord( irecord, orecord, weigth );
	analyse_Muon		->analyseRecord( irecord, orecord, weigth );
	analyse_Electron	->analyseRecord( irecord, orecord, weigth );
	analyse_Photon		->analyseRecord( irecord, orecord, weigth );
	analyse_Jet		->analyseRecord( irecord, orecord, weigth );
	analyse_Mis		->analyseRecord( irecord, orecord, weigth );
	analyse_BJet		->analyseRecord( irecord, orecord, weigth );
	analyse_CJet		->analyseRecord( irecord, orecord, weigth );
	analyse_Tau		->analyseRecord( irecord, orecord, weigth );
	analyse_Calibration	->analyseRecord( irecord, orecord, weigth );
}

/*
 * Print info about AcerDET class
 */
void AcerDET::printInfo() const {
	// AcerDET version number and release date
	printf ("**********************************************************\n");
	printf ("*                  AcerDET, version: 2.0                 *\n");
	printf ("*                 Released at: 30.06.2015                *\n");
	printf ("*                                                        *\n");
	printf ("*  Simplied event simulation and reconstruction package  *\n");
	printf ("*                                                        *\n");
	printf ("*           by E. Richter-Was (Institute of Physics)     *\n");
	printf ("*          and P. Mikos (Theoretical Computer Science)   *\n");
	printf ("*         Jagiellonian University, Cracow, Poland        *\n");
	printf ("**********************************************************\n");
	printf ("\n");
	
	// info about subclasses
	printf (" Initial configuration:\n");
	analyse_Cell		->printInfo();
	analyse_Cluster		->printInfo();
	analyse_Muon		->printInfo();
	analyse_Electron	->printInfo();
	analyse_Photon		->printInfo();
	analyse_Jet		->printInfo();
	analyse_Mis		->printInfo();
	analyse_BJet		->printInfo();
	analyse_CJet		->printInfo();
	analyse_Tau		->printInfo();
	analyse_Calibration	->printInfo();
}

void AcerDET::printResults() const {
	analyse_Cell		->printResults();
	analyse_Cluster		->printResults();
	analyse_Muon		->printResults();
	analyse_Electron	->printResults();
	analyse_Photon		->printResults();
	analyse_Jet		->printResults();
	analyse_Mis		->printResults();
	analyse_BJet		->printResults();
	analyse_CJet		->printResults();
	analyse_Tau	       	->printResults();
	analyse_Calibration	->printResults();
}

void AcerDET::storeHistograms() {
	histos->storeHistograms();
}
