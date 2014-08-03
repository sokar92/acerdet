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
	
	histo_cJets		("CJet: c-jets multiplicity", 0.0, 10.0, 10),
	histo_cQuarks	("CJet: c-quarks HARD multiplicity", 0.0, 10.0, 10),
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
		
		if (part.type == PT_CJET) {  // && part.statusID != 21 ? po co
			// if there is a c-quark found before hadronization
			// if there are still jets
			if (!orecord.Jets.isEmpty()) {
				Bool_t CJET = true;
				Int32_t JETC = 0;
				
				// and this c-quark is the last one in the FSR cascade
				if (part.hasDaughter()) {
					for (int j=part.daughters.first; j<=part.daughters.second; ++j) {
						if (parts[j].type == PT_CJET) 
							CJET = false;
					}
				}

				if (!CJET) 
					continue;

				PT = part.pT();
				if (PT < PTCMIN) 
					continue;

				ETA = part.getEta(); 
				if (abs(ETA) > ETCMAX)
					continue;

				PHI = part.getPhi();
				
				// mark c-jet
				DR = 100.0;
				for (int j=0; j<orecord.Jets.size(); ++j) {
					if (orecord.Jets[j].type != PT_BJET) {
						DDR = sqrt(
							pow(ETA - orecord.Jets[j].eta_rec (PJET(II,3)), 2) +
							pow(PHI - orecord.Jets[j].phi_rec (PJET(II,4)), 2)
						);
						
						if (abs(PHI - orecord.Jets[j].phi_rec (PJET(II,4)) ) > PI)
							DDR = sqrt(
								pow(ETA - PJET(II,3), 2) +
								pow(abs(PHI - PJET(II,4))-2*PI, 2)
							);

						if (DDR < DR) 
							JETC = j;
						
						DR = min(DDR,DR);
					}
				}

				if (DR > RJC) {
					continue;
				}

				// labell  c-jet
				KJET(JETC,2) = 4							// JETC = indeks znalezionego jetu
				KJET(JETC,5) = I
				NJETC = NJETC + 1							// kolejny cjet znaleziony
			}
		}
	}

	// store count in histogram
	histo_cJets.insert(NJETC);
	
	// check partons
	Int32_t IQUAC = 0, ICJET = 0;
	for (int i=6; i<=NSTOP; ++i) {
		const Particle& part = parts[i];
		
		if (part.type == PT_CJET) {
			PT = part.pT();
			ETA = part.getEta(); 
			PHI = part.getPhi();

			if (abs(ETA) < ETCMAX && PT > ETJET) {
				IQUAC++;
				DR = 18.0;

				for (int j=0; j<orecord.Jets.size(); ++j) {
					if (orecord.Jets[j].type == PT_CJET) {  // skoro iterujemy sie po CJetach to po co sprawdzac to jeszcze raz
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
							ICJET = j;
							
						if (DDR < DR) 
							DR = DDR;
					}
				}
			}
		}
	}

	// store quarks-count in histogram
	histo_cQuarks.insert(IQUAC);

	for (int i=0; i<orecord.Jets.size(); ++i) {
		PTREC = 0;
		DETRMIN = RCONE;

		if (orecord.Jets[i].type == PT_CJET) {
			for (int j=6; j<parts.size(); ++j) {
				const Particle& part = parts[j];
				
				if (part.stateID != 21 || part.type != PT_CJET) 
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
				
				if (DETR > DETRMIN) 
					continue;

				PTREC = PT;
				DETRMIN = DETR;
			}
		}

		if (PTREC) {
			histo_delta.insert(DETRMIN);
			histo_pT.insert(PJET(IJET,5) / PTREC);
		}
	}
	*/
}

void CJet::printResults() const {
	if (KEYBCL) {
		printf ("***********************************\n");
		printf ("*                                 *\n");
		printf ("*     ***********************     *\n");
		printf ("*     ***   Output from   ***     *\n");
		printf ("*     ***  analyse::CJet  ***     *\n");
		printf ("*     ***********************     *\n");
		printf ("*                                 *\n");
		printf ("***********************************\n");
	
		printf (" Analysed records: %d\n", IEVENT);
		histo_cJets		.print( true ); // IDENT + 11
		histo_cQuarks	.print( true ); // IDENT + 21
		histo_delta		.print( true ); // IDENT + 23
		histo_pT		.print( true ); // IDENT + 24
	}
}
