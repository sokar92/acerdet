#include "Muon.h"
#include <cstdio>

using namespace AcerDet::analyse;

Muon::Muon( const Configuration& config ) :
	ETCLU	( config.Cluster.MinEt ),
	RCONE	( config.Cluster.ConeR ),

	PTMUMIN	( config.Muon.MinMomenta ),
	ETAMAX	( config.Muon.MaxEta ),
	RISOLJ	( config.Muon.MinIsolRlj ),
	RDEP	( config.Muon.ConeR ),
	EDMAX	( config.Muon.MaxEnergy ),
	PTMUMINT( config.Muon.MinMomenta ),    // PROBLEM!!!

	KEYHID	( config.Flag.HistogramId ),
	KEYSME	( config.Flag.Smearing ),

	IEVENT	( 0 ),
	
	histo_nonisol	("Muon: non-isolated", 0.0, 10.0, 10),
	histo_isol		("Muon: isolated", 0.0, 10.0, 10),
	histo_hard		("Muon: hard", 0.0, 10.0, 10),
	histo_sum		("Muon: hard+isol", 0.0, 10.0, 10)
{}

Muon::~Muon() {}

void Muon::printInfo() const {
	// print out title
	printf ("*************************************\n");
	printf ("*                                   *\n");
	printf ("*     *************************     *\n");
	printf ("*     ***   analyse::Muon   ***     *\n");
	printf ("*     *************************     *\n");
	printf ("*                                   *\n");
	printf ("*************************************\n");

	// print out basic params
	printf ("\n\t... muon isolation ...\n");
	printf ("min. muon p_T %lf\n", PTMUMIN);
	printf ("max. muon eta %lf\n", ETAMAX);
	printf ("min R_lj for isolation %lf\n", RISOLJ);
	printf ("R for energy deposit %lf\n", RDEP);
	printf ("max E_dep for isolation %lf\n", EDMAX);
	printf ("min muon p_T unsmea %lf\n", PTMUMINT);
	printf ("smearing %s\n", KEYSME ? "on" : "off");
	printf ("\n");
}

void Muon::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	//printf ("Muon: analyse record\n");
}
