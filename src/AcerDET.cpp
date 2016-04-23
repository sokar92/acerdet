#include "AcerDET.h"
#include <cstdio>

#include "analyse/BJet.h"
#include "analyse/Calibration.h"
#include "analyse/Cell.h"
#include "analyse/CJet.h"
#include "analyse/Cluster.h"
#include "analyse/Electron.h"
#include "analyse/Jet.h"
#include "analyse/Mis.h"
#include "analyse/Muon.h"
#include "analyse/Photon.h"
#include "analyse/Tau.h"
#include "analyse/Test_Histograms.h"

#include "analyse/FastJet_Clustering.h"

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
	
	analyse_phases.push_back( new analyse::Cell( config, histoManager, partProvider ) );
	
	if (config.Flag.UseFastJet)
		analyse_phases.push_back( new analyse::FastJet_Clustering( config, histoManager, partProvider ) );
	else
		analyse_phases.push_back( new analyse::Cluster( config, histoManager, partProvider ) );
	
	analyse_phases.push_back( new analyse::Muon( config, histoManager ) );
	analyse_phases.push_back( new analyse::Electron( config, histoManager ) );
	analyse_phases.push_back( new analyse::Photon( config, histoManager ) );
	analyse_phases.push_back( new analyse::Jet( config, histoManager, partProvider ) );
	analyse_phases.push_back( new analyse::Mis( config, histoManager ) );
	analyse_phases.push_back( new analyse::BJet( config, histoManager ) );
	analyse_phases.push_back( new analyse::CJet( config, histoManager ) );
	analyse_phases.push_back( new analyse::Tau( config, histoManager ) );
	analyse_phases.push_back( new analyse::Calibration( config, histoManager ) );
	analyse_phases.push_back( new analyse::Test_Histograms( config, histoManager ) );
}

/*
 * Destructor - not used
 */
AcerDET::~AcerDET() {
	histos = NULL;
	
	for (int i=0; i<analyse_phases.size(); i++)
		delete analyse_phases[i];
}

/*
 * Analyse single InputRecord from event
 */
void AcerDET::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord, Real64_t weigth ) {

	if (!histos_initialized) {
		histos->init();
		histos_initialized = true;
	}
	
	for (std::vector<analyse::IAnalysePhase*>::const_iterator it = analyse_phases.begin(); it != analyse_phases.end(); it++)
		(*it)->analyseRecord( irecord, orecord, weigth );
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
	
	for (std::vector<analyse::IAnalysePhase*>::const_iterator it = analyse_phases.begin(); it != analyse_phases.end(); it++)
		(*it)->printInfo();
}

void AcerDET::printResults() const {
	for (std::vector<analyse::IAnalysePhase*>::const_iterator it = analyse_phases.begin(); it != analyse_phases.end(); it++)
		(*it)->printResults();
}

void AcerDET::storeHistograms() {
	histos->storeHistograms();
}
