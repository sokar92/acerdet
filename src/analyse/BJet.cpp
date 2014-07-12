#include "BJet.h"
#include <cstdio>

using namespace AcerDet::analyse;

BJet::BJet( const Configuration& config ) :
	ETJET	( config.Jet.MinEnergy ),
	RCONE	( config.Cluster.ConeR ),

	PTBMIN	( config.BJet.MinMomenta ),
	ETBMAX	( config.BJet.MaxEta ),
	RJB	( config.BJet.MaxRbj ),

	KEYHID	( config.Flag.HistogramId ),
	KEYBCL	( config.Flag.BCJetsLabeling ),

	IEVENT	( 0 ),
	
	histo_bJets		("BJet: b-jets multiplicity", 0.0, 10.0, 10),
	histo_bQuarks	("BJet: b-quarks HARD multiplicity", 0.0, 10.0, 10),
	histo_delta		("BJet: delta r bjet-bquark", 0.0, 0.5, 50),
	histo_pT		("BJet: pTbjet / pTbquark", 0.0, 2.0, 50)
{}

BJet::~BJet() {}

void BJet::printInfo() const {
	// print out title
	printf ("*************************************\n");
	printf ("*                                   *\n");
	printf ("*     *************************     *\n");
	printf ("*     ***   analyse::BJet   ***     *\n");
	printf ("*     *************************     *\n");
	printf ("*                                   *\n");
	printf ("*************************************\n");

	// print out basic params
	printf ("\n\t... jets labeling ...\n");
	printf ("labeling on/off %s\n", KEYBCL ? "on" : "off");
	if (KEYBCL) {
		printf ("\tbjets ...\n");
		printf ("min b-quark p_T %lf\n", PTBMIN);
		printf ("max b-quark eta %lf\n", ETBMAX);
		printf ("max R_bj for b-jets %lf\n", RJB);
	}
	printf ("\n");
}

#define KEY_HISTO 21

void BJet::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	//printf ("BJet: analyse record\n");
	
	// option is off
	if(!KEYBCL)
		return;
		
	// next analysed event
	IEVENT++;
	
	// search for history partition
	const vector<Particle>& parts = irecord.particles();
	int n = parts.size();
	
	int nstart = 0;
	while(nstart < n && parts[nstart].state != PS_HISTORY) 
		nstart++;
	int nstop = nstart-1;
}
