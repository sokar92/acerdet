#include "Calibration.h"
#include <cstdio>

using namespace AcerDet::analyse;

Calibration::Calibration( const Configuration& config ) :
	RCONE	( config.Cluster.ConeR ),
	
	KEYHID	( config.Flag.HistogramId ),
	KEYCAL	( config.Flag.JetCalibration ),

	IEVENT	( 0 ),
	IDENT	( 1000 + KEYHID )
{
}

Calibration::~Calibration() {
}

void Calibration::printInfo() const {
	// print out title
	printf ("********************************************\n");
	printf ("*                                          *\n");
	printf ("*     ********************************     *\n");
	printf ("*     ***   analyse::Calibration   ***     *\n");
	printf ("*     ********************************     *\n");
	printf ("*                                          *\n");
	printf ("********************************************\n");

	// print out basic params
	printf ("\n\t jets calibration ....\n");
	printf (" calibration %s\n", KEYCAL ? "on" : "off");
	printf ("\n");
}

void Calibration::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	//printf ("Calibration: analyse record\n");
}
