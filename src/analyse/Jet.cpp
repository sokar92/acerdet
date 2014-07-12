#include "Jet.h"
#include <cstdio>

using namespace AcerDet::analyse;

Jet::Jet( const Configuration& config ) :
	ETJET	( config.Jet.MinEnergy ),
	ETAJET	( config.Jet.RapidityCoverage ),
	RCONE	( config.Cluster.ConeR ),
	PTMIN	( config.Cell.MinpT ),
	CALOTH	( config.Cell.EtaTransition ),

	KEYHID	( config.Flag.HistogramId ),
	KEYSME	( config.Flag.Smearing ),
	KEYFLD	( config.Flag.BField ),
	KFINVS	( config.Flag.SusyParticle ),

	IEVENT	( 0 ),
	
	histo_bJets				("Jet: multiplicity", 0.0, 10.0, 10),
	histo_delta_phi			("Jet: delta phi jet-barycentre", -0.5, 0.5, 50),
	histo_delta_eta			("Jet: delta eta jet-barycentre", -0.5, 0.5, 50),
	histo_delta_barycenter	("Jet: delta r jet-barycentre", 0.0, 0.5, 50),
	histo_delta_parton		("Jet: delta r jet-parton", 0.0, 0.5, 50),
	histo_pT_bySum			("Jet: pTjet / SumpTParticle", 0.0, 2.0, 50),
	histo_pT_byPart			("Jet: pTjet / pTparton", 0.0, 2.0, 50)
{}

Jet::~Jet() {}

void Jet::printInfo() const {
	// print out title
	printf ("************************************\n");
	printf ("*                                  *\n");
	printf ("*     ************************     *\n");
	printf ("*     ***   analyse::Jet   ***     *\n");
	printf ("*     ************************     *\n");
	printf ("*                                  *\n");
	printf ("************************************\n");

	// print out basic params
	printf ("\n\t clusters definition ....\n");
	printf (" R cone %lf\n", RCONE);
	printf ("\t jets definition ....\n");
	printf (" E_T_jets [GeV] %lf\n", ETJET);
	printf (" eta coverage jets %lf\n", ETAJET);
    printf (" smearing %s\n", KEYSME ? "on" : "off");
    printf (" B-field %s\n", KEYFLD ? "on" : "off");
	printf ("\n");
}

void Jet::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	//printf ("Jet: analyse record\n");
}
