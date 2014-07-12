#include "Photon.h"
#include <cstdio>

using namespace AcerDet::analyse;

Photon::Photon( const Configuration& config ) :
	ETCLU	( config.Cluster.MinEt ),
	RCONE	( config.Cluster.ConeR ),

	PTLMIN	( config.Photon.MinMomenta ),
	ETAMAX	( config.Photon.MaxEta ),
	RJE	( config.Photon.MinJetsRlj ),
	RISOLJ	( config.Photon.MinIsolRlj ),
	RDEP	( config.Photon.ConeR ),
	EDMAX	( config.Photon.MaxEnergy ),

	KEYHID	( config.Flag.HistogramId ),
	KEYSME	( config.Flag.Smearing ),

	IEVENT	( 0 ),
	
	histo_isol	("Photon: multiplicity isolated", 0.0, 10.0, 10),
	histo_hard	("Photon: multiplicity hard", 0.0, 10.0, 10),
	histo_sum	("Photon: multipliciyt isol+hard", 0.0, 10.0, 10)
{}

Photon::~Photon() {}

void Photon::printInfo() const {
	// print out title
	printf ("***************************************\n");
	printf ("*                                     *\n");
	printf ("*     ***************************     *\n");
	printf ("*     ***   analyse::Photon   ***     *\n");
	printf ("*     ***************************     *\n");
	printf ("*                                     *\n");
	printf ("***************************************\n");

	// print out basic params
	printf ("\n\t... photon isolation ...\n");
	printf ("min. photon p_T %lf\n", PTLMIN);
	printf ("max. photon eta %lf\n", ETAMAX);
	printf ("max R_gam-clust %lf\n", RJE);
	printf ("min R_isol %lf\n", RISOLJ);
	printf ("R for energy deposit %lf\n", RDEP);
	printf ("max E_dep for isolation %lf\n", EDMAX);
	printf ("smearing %s\n", KEYSME ? "on" : "off");
	printf ("\n");

}

void Photon::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	//printf ("Photon: analyse record\n");
}
