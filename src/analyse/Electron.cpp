#include "Electron.h"
#include <cstdio>

using namespace AcerDet::analyse;

Electron::Electron( const Configuration& config ) :
	ETCLU	( config.Cluster.MinEt ),
	RCONE	( config.Cluster.ConeR ),

	PTLMIN	( config.Electron.MinMomenta ),
	ETAMAX	( config.Electron.MaxEta ),
	RJE	( config.Electron.MinJetsRlj ),
	RISOLJ	( config.Electron.MinIsolRlj ),
	RDEP	( config.Electron.ConeR ),
	EDMAX	( config.Electron.MaxEnergy ),

	KEYHID	( config.Flag.HistogramId ),
	KEYSME	( config.Flag.Smearing ),

	IEVENT	( 0 ),
	
	histo_isol	("Electron: multiplicity isolated", 0.0, 10.0, 10),
	histo_hard	("Electron: multiplicity hard", 0.0, 10.0, 10),
	histo_sum	("Electron: multipliciyt isol+hard", 0.0, 10.0, 10)
{}

Electron::~Electron() {}

void Electron::printInfo() const {
	// print out title
	printf ("*****************************************\n");
	printf ("*                                       *\n");
	printf ("*     *****************************     *\n");
	printf ("*     ***   analyse::Electron   ***     *\n");
	printf ("*     *****************************     *\n");
	printf ("*                                       *\n");
	printf ("*****************************************\n");

	// print out basic params
	printf ("\n\t... electron isolation ...\n");
	printf ("min. lepton p_T %lf\n", PTLMIN);
	printf ("max. lepton eta %lf\n", ETAMAX);
	printf ("max R_ej for ele-clu %lf\n", RJE);
	printf ("min R_lj for isolation %lf\n", RISOLJ);
	printf ("R for energy deposit %lf\n", RDEP);
	printf ("max E_dep for isolation %lf\n", EDMAX);
	printf ("smearing %s\n", KEYSME ? "on" : "off");
	printf ("\n");

}

void Electron::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	//printf ("Electron: analyse record\n");
}
