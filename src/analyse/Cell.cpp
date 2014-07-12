#include "Cell.h"
#include <cstdio>

using namespace AcerDet::analyse;

Cell::Cell( const Configuration& config ) :
	ETACEL	( config.Cell.RapidityCoverage ),
	PTMIN	( config.Cell.MinpT ),
	ETTHR	( config.Cell.MinEt ),
	CALOTH	( config.Cell.EtaTransition ),
	DBETA	( config.Cell.GranularityEta ),
	DBPHI	( config.Cell.GranularityPhi ),
	
	KEYHID	( config.Flag.HistogramId ),
	KEYFLD	( config.Flag.BField ),
	KFINVS	( config.Flag.SusyParticle ),

	IEVENT	( 0 ),
	
	histo ("Cell: multiplicity", 0.0, 500.0, 50)
{}

Cell::~Cell() {}

void Cell::printInfo() const {
	// print out title
	printf ("*************************************\n");
	printf ("*                                   *\n");
	printf ("*     *************************     *\n");
	printf ("*     ***   analyse::Cell   ***     *\n");
	printf ("*     *************************     *\n");
	printf ("*                                   *\n");
	printf ("*************************************\n");

	// print out basic params
	printf ("\n\t clusters definition ...\n");
	printf (" eta coverage  %lf\n", ETACEL);
	printf (" E_T_min cell thresh %lf\n", ETTHR);
	printf (" eta gran. transition %lf\n", CALOTH);
	printf (" gran in eta(central) %lf\n", DBETA);
	printf (" gran in phi(central) %lf\n", DBPHI);
	printf ("\n\t B field apply ....\n");
	printf (" B-field %s\n", KEYFLD ? "on" : "off");
	if (KEYFLD) {
		printf (" p_T min non looping %lf\n", PTMIN);
	}
	printf ("\n");
}

void Cell::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	//printf ("Cell: analyse record\n");
}
