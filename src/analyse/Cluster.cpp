#include "Cluster.h"
#include <cstdio>

using namespace AcerDet::analyse;

Cluster::Cluster( const Configuration& config ) :
	ETCLU	( config.Cluster.MinEt ),
	RCONE	( config.Cluster.ConeR ),
	ETACLU	( config.Cluster.RapidityCoverage ),
	ETINI	( config.Cluster.MinEtInit ),
	
	PTMIN	( config.Cell.MinpT ),
	CALOTH	( config.Cell.EtaTransition ),
	
	KEYHID	( config.Flag.HistogramId ),
	KEYFLD	( config.Flag.BField ),
	KFINVS	( config.Flag.SusyParticle ),

	IEVENT	( 0 ),
	
	histo_bJets				("Cluster: multiplicity", 0.0, 10.0, 10),
	histo_delta_phi			("Cluster: delta phi clu-barycentre", -0.5, 0.5, 50),
	histo_delta_eta			("Cluster: delta eta clu-barycentre", -0.5, 0.5, 50),
	histo_delta_barycenter	("Cluster: delta r clu-barycentre", 0.0, 0.5, 50),
	histo_delta_parton		("Cluster: delta r clu-parton", 0.0, 0.5, 50),
	histo_pT_bySum			("Cluster: pTclu / SumpTParticle", 0.0, 2.0, 50),
	histo_pT_byPart			("Cluster: pTclu / pTparton", 0.0, 2.0, 50)
{}

Cluster::~Cluster() {}

void Cluster::printInfo() const {
	// print out title
	printf ("****************************************\n");
	printf ("*                                      *\n");
	printf ("*     ****************************     *\n");
	printf ("*     ***   analyse::Cluster   ***     *\n");
	printf ("*     ****************************     *\n");
	printf ("*                                      *\n");
	printf ("****************************************\n");

	// print out basic params
	printf ("\n\t clusters definition ...\n");
	printf (" E_T_min cluster %lf\n", ETCLU);
	printf (" E_T_min cell initia %lf\n", ETINI);
	printf (" R cone %lf\n", RCONE);
	printf (" eta coverage %lf\n", ETACLU);
	printf (" eta gran. transition %lf\n", CALOTH);
	printf ("\n");
}

void Cluster::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	//printf ("Cluster: analyse record\n");
}
