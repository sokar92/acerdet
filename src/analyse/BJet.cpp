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

void BJet::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	if (!KEYBCL)
		return;
/*
	// new event to compute
	IEVENT++;

	// reference to particles container
	const vector<Particle>& parts = irecord.particles();
	Int32_t N = parts.size();

	// znajdz poczatek danych
	Int32_t NSTOP = 0, NSTART = 1;
	for (int i=0; i<parts.size(); ++i) {
		if (parts[i].stateID != 21) {
			NSTOP = i-1;
			NSTART = i;
			break;
		}
	}
	
	// look for c-jets
	for (int i=NSTART; i<parts.size(); ++i) {
		const Particle& part = parts[i];
		
		if (part.type == PT_BJET) {  // && part.statusID != 21 ? po co
			// if there is a b-quark found before hadronization
			// if there are still jets
			if (!orecord.Jets.isEmpty()) {
				Bool_t BJET = true;
				Int32_t JETB = -1;
				
				// and this b-quark is the last one in the FSR cascade
				if (part.hasDaughter()) {
					for (int j=part.daughters.first; j<=part.daughters.second; ++j) {
						if (parts[j].type == PT_BJET)
							BJET = false;
					}
				}

				if (!BJET) 
					continue;

				PT = part.pT(); 
				if (PT < PTBMIN) 
					continue;

				ETA = part.getEta();
				if (abs(ETA) > ETBMAX)
					continue;
			
				PHI = part.getPhi();
				
				// mark b-jet
				DR = 100.0;
				for (int j=0; j<orecord.Jets.size(); ++j) {
					DDR = sqrt(
						pow(ETA - PJET(II,3), 2) +
						pow(PHI - PJET(II,4), 2)
					);
					
					if (abs(PHI - PJET(II,4)) > PI )
						DDR = sqrt(
							pow(ETA - PJET(II,3), 2) + 
							pow(abs(PHI - PJET(II,4))-2*PI, 2)
						);
						
					if (DDR < DR) 
						JETB = j;
						
					DR = min(DDR,DR);
				}

				if (DR > RJB)
					continue;
					
				// label  b-jet
				KJET(JETB,2) = 5;
				KJET(JETB,5) = I;
				NJETB = NJETB + 1;
			}
		}
	}

	histo_bJets.insert(NJETB);
	
	// check partons
	Int32_t IQUAB = 0, IBJET = 0;
	for (int i=6; i<=NSTOP; ++i) {								// 7? -> barcode = numery czastek w zdarzeniu
		const Particle& part = parts[i];
		
		if (part.type == PT_BJET) {
			PT = part.pT();
			ETA = part.getEta();
			PHI = part.getPhi();

			if (abs(ETA) < ETBMAX && PT > ETJET) {
				IQUAB++;
				DR = 18.0;

				for (int j=0; j<orecord.Jets.size(); ++j) {
					if (orecord.Jets[j].type == PT_BJET) {
						DDR = sqrt(
							pow(ETA - PJET(II,3), 2) +
							pow(PHI - PJET(II,4), 2)
						);
						
						if (abs(PHI - PJET(II,4)) > PI)
							DDR = sqrt(
								pow(ETA - PJET(II,3), 2) +
								pow(abs(PHI - PJET(II,4))-2*PI, 2)
							);
						
						if (DDR < DR) 
							IBJET = j;
							
						if (DDR < DR) 
							DR = DDR;
					}
				}
			}
		}
	}

	// store count of b-quarks in histogram
	histo_bQuarks.insert(IQUAB);
	
	for (int i=0; i<orecord.Jets.size(); ++i) {
		PTREC = 0;
		DETRMIN = RCONE;

		if (orecord.Jets[i].type == 5) {
			for (int j=6; j<parts.size(); ++j) {
				const Particle& part = parts[j];
				
				if (part.type != PT_BJET) // && part.stateID != 21
					continue;

				PT = part.pT();
				ETA = part.getEta();
				PHI = part.getPhi();
				DPHIA = abs(orecord.Jets[i].phi (PJET(IJET,4)) - PHI); 

				if (DPHIA > PI) 
					DPHIA -= 2*PI;

				DETR = sqrt(
					pow(ETA - orecord.Jets[i].eta (PJET(IJET,3)) , 2) +
					pow(DPHIA, 2)
				);
				
				if(DETR < DETRMIN) {
					PTREC = PT;
					DETRMIN = DETR;
				}
			}
		}

		if (PTREC != 0) {
			histo_delta.insert(DETRMIN);
			histo_pT.insert(orecord.Jets[i].pT (PJET(IJET,5)) / PTREC);
		}
	}
	*/
}

void BJet::printResults() const {
	if (KEYBCL) {
		printf ("***********************************\n");
		printf ("*                                 *\n");
		printf ("*     ***********************     *\n");
		printf ("*     ***   Output from   ***     *\n");
		printf ("*     ***  analyse::BJet  ***     *\n");
		printf ("*     ***********************     *\n");
		printf ("*                                 *\n");
		printf ("***********************************\n");
	
		printf (" Analysed records: %d\n", IEVENT);
		histo_bJets		.print( true );
		histo_bQuarks	.print( true );
		histo_delta		.print( true );
		histo_pT		.print( true );
	}
}
