#include "CJet.h"
#include <cstdio>

using namespace AcerDet::analyse;

CJet::CJet( const Configuration& config ) :
	ETJET	( config.Jet.MinEnergy ),
	RCONE	( config.Cluster.ConeR ),

	PTCMIN	( config.CJet.MinMomenta ),
	ETCMAX	( config.CJet.MaxEta ),
	RJC	( config.CJet.MaxRcj ),

	KEYHID	( config.Flag.HistogramId ),
	KEYBCL	( config.Flag.BCJetsLabeling ),

	IEVENT	( 0 ),
	
	histo_bJets		("CJet: c-jets multiplicity", 0.0, 10.0, 10),
	histo_bQuarks	("CJet: c-quarks HARD multiplicity", 0.0, 10.0, 10),
	histo_delta		("CJet: delta r cjet-cquark", 0.0, 0.5, 50),
	histo_pT		("CJet: pTcjet / pTcquark", 0.0, 2.0, 50)
{}

CJet::~CJet() {}

void CJet::printInfo() const {
	// print out title
	printf ("*************************************\n");
	printf ("*                                   *\n");
	printf ("*     *************************     *\n");
	printf ("*     ***   analyse::CJet   ***     *\n");
	printf ("*     *************************     *\n");
	printf ("*                                   *\n");
	printf ("*************************************\n");

	// print out basic params
	printf ("\n\t... jets labeling ...\n");
	printf ("labeling on/off %s\n", KEYBCL ? "on" : "off");
	if (KEYBCL) {
		printf ("\tcjets ...\n");
		printf ("min c-quark p_T %lf\n", PTCMIN);
		printf ("max c-quark eta %lf\n", ETCMAX);
		printf ("max R_cj for c-jets %lf\n", RJC);
	}
	printf ("\n");
}

void CJet::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	//printf ("CJet: analyse record\n");
}
