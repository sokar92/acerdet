#include "Mis.h"
#include <cstdio>

using namespace AcerDet::analyse;

Mis::Mis( const Configuration& config ) :
	PTMUMIN	( config.Muon.MinMomenta ),
	ETAMAX	( config.Muon.MaxEta ),
	ETCELL	( config.Tau.MinpT ),
	CALOTH	( config.Cell.EtaTransition ),
	
	KEYHID	( config.Flag.HistogramId ),
	KEYSME	( config.Flag.Smearing ),
	KFINVS	( config.Flag.SusyParticle ),

	IEVENT	( 0 ),

	histo_reconstructed_pT			("Mis: reconstructed p_T", 0.0, 200.0, 50),
	histo_reconstructed_pT_cells	("Mis: reconstructed p_T + cells", 0.0, 200.0, 50),
	histo_pTmiss					("Mis: pTmiss", 0.0, 200.0, 50),
	histo_pTnu						("Mis: pT nu", 0.0, 200.0, 50)
{}

Mis::~Mis() {}

void Mis::printInfo() const {
	// print out title
	printf ("************************************\n");
	printf ("*                                  *\n");
	printf ("*     ************************     *\n");
	printf ("*     ***   analyse::Mis   ***     *\n");
	printf ("*     ************************     *\n");
	printf ("*                                  *\n");
	printf ("************************************\n");

	// print out basic params
	printf ("\n\t muon coverage\n");
	printf (" min. muon p_T %lf\n", PTMUMIN);
	printf (" max. muon eta %lf\n", ETAMAX);
	printf ("\t unused cells ...\n");
	printf (" smearing %s\n", KEYSME ? "on" : "off");
	printf (" cells threshold %lf\n", ETCELL);
	printf ("\t invisible particles ...\n");
	printf (" KF code for invis %d\n", KFINVS);
	printf ("\n");
}

void Mis::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	//printf ("Mis: analyse record\n");
}

void Mis::printResults() const {
	printf ("***********************************\n");
	printf ("*                                 *\n");
	printf ("*     ***********************     *\n");
	printf ("*     ***   Output from   ***     *\n");
	printf ("*     ***  analyse::Miss  ***     *\n");
	printf ("*     ***********************     *\n");
	printf ("*                                 *\n");
	printf ("***********************************\n");
	
	printf (" Analysed records: %d\n", IEVENT);
	histo_reconstructed_pT			.print( true );
	histo_reconstructed_pT_cells	.print( true );
	histo_pTmiss					.print( true );
	histo_pTnu						.print( true );
}
