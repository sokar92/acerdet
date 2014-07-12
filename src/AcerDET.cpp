#include "AcerDET.h"
#include <cstdio>

using namespace AcerDet;

/*
 * Construct AcerDET object from given configuration
 */
AcerDET::AcerDET( const conf::Configuration& config ) : 
		analyse_BJet		( config ),
		analyse_Calibration	( config ),
		analyse_Cell		( config ),
		analyse_CJet		( config ),
		analyse_Cluster		( config ),
		analyse_Electron	( config ),
		analyse_Jet			( config ),
		analyse_Mis			( config ),
		analyse_Muon		( config ),
		analyse_Photon		( config ),
		analyse_Tau			( config )
{
}

/*
 * Destructor - not used
 */
AcerDET::~AcerDET() {
}

/*
 * Analyse single InputRecord from event (dodatkowy parametr pointer na kontener)
 */
void AcerDET::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	analyse_Cell		.analyseRecord( irecord, orecord );
	analyse_Cluster		.analyseRecord( irecord, orecord );
	analyse_Muon		.analyseRecord( irecord, orecord );
	analyse_Electron	.analyseRecord( irecord, orecord );
	analyse_Photon		.analyseRecord( irecord, orecord );
	analyse_Jet			.analyseRecord( irecord, orecord );
	analyse_Mis			.analyseRecord( irecord, orecord );
	analyse_BJet		.analyseRecord( irecord, orecord );
	analyse_CJet		.analyseRecord( irecord, orecord );
	analyse_Tau			.analyseRecord( irecord, orecord );
	analyse_Calibration	.analyseRecord( irecord, orecord );
}

/*
 * Print info about AcerDET class
 */
void AcerDET::printInfo() const {
	// AcerDET version number and release date
	printf ("**********************************************************\n");
	printf ("*                  AcerDET, version: 2.0                 *\n");
	printf ("*                 Released at: xx.yy.2014                *\n");
	printf ("*                                                        *\n");
	printf ("*  Simplied event simulation and reconstruction package  *\n");
	printf ("*                                                        *\n");
	printf ("*              by E. Richter-Was & P. Mikos              *\n");
	printf ("*             Institute of Computer Science              *\n");
	printf ("*         Jagiellonian University, Cracow, Poland        *\n");
	printf ("**********************************************************\n");
	printf ("\n");
	
	// info about subclasses
	printf (" Initial configuration:\n");
	analyse_Cell		.printInfo();
	analyse_Cluster		.printInfo();
	analyse_Muon		.printInfo();
	analyse_Electron	.printInfo();
	analyse_Photon		.printInfo();
	analyse_Jet			.printInfo();
	analyse_Mis			.printInfo();
	analyse_BJet		.printInfo();
	analyse_CJet		.printInfo();
	analyse_Tau			.printInfo();
	analyse_Calibration	.printInfo();
}
