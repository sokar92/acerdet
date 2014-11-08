#include "CJet.h"
#include <cstdio>
using namespace AcerDet::analyse;

#include "../core/Functions.h"
using namespace AcerDet::core;

CJet::CJet( const Configuration& config, IHistogramManager *histoMng ) :
	ETJET	( config.Jet.MinEnergy ),
	RCONE	( config.Cluster.ConeR ),

	PTCMIN	( config.CJet.MinMomenta ),
	ETCMAX	( config.CJet.MaxEta ),
	RJC	( config.CJet.MaxRcj ),

	KEYHID	( config.Flag.HistogramId ),
	KEYBCL	( config.Flag.BCJetsLabeling ),

	IEVENT	( 0 ),

	histoManager(histoMng),
	histoRegistered( false )
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
	
	Int32_t idhist = 800 + KEYHID;
	if (!histoRegistered) {
		histoRegistered = true;
		histoManager
			->registerHistogram(idhist+11, "CJet: c-jets multiplicity", 10, 0.0, 10.0);
		histoManager
			->registerHistogram(idhist+21, "CJet: c-quarks HARD multiplicity", 10, 0.0, 10.0);
		histoManager
			->registerHistogram(idhist+23, "CJet: delta r cjet-cquark", 50, 0.0,  5.0);
		histoManager
			->registerHistogram(idhist+24, "CJet: pTcjet/pTcquark", 50, 0.0,  2.0);
	}

	// do not use this algorithm
	if (!KEYBCL)
		return;

	// variables
	Real64_t PT, PHI, ETA, DR, DDR, DETR, DPHIA, DETRMIN, PTREC;

	// new event to compute
	IEVENT++;

	// reference to particles container
	const vector<Particle>& parts = irecord.particles();
	
	// look for c-jets
	for (int i=0; i<parts.size(); ++i) {
		const Particle& part = parts[i];
		
		if (part.status == PS_FINAL 
		&& part.type == PT_CJET) {
			// if there is a c-quark found before hadronization
			// if there are still jets
			if (!orecord.Jets.empty()) {
				Bool_t CJET = true;
				Int32_t JETC = -1;
				
				// and this c-quark is the last one in the FSR cascade
				if (part.hasDaughter()) {
					for (int j=part.daughters.first; j<=part.daughters.second; ++j) {
						if (parts[j-1].type == PT_CJET) 
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
							pow(ETA - orecord.Jets[j].eta_rec, 2) +
							pow(PHI - orecord.Jets[j].phi_rec, 2)
						);
						
						if (abs(PHI - orecord.Jets[j].phi_rec) > PI)
							DDR = sqrt(
								pow(ETA - orecord.Jets[j].eta_rec, 2) +
								pow(abs(PHI - orecord.Jets[j].phi_rec) - 2*PI, 2)
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
				// KJET(JETC,2) = 4							// JETC = indeks znalezionego jetu
				// KJET(JETC,5) = I
				// NJETC = NJETC + 1						// kolejny cjet znaleziony
			}
		}
	}

	// store count in histogram TODO (replace NJETC when cJets was in output record)
	// histoManager
	//   ->insert(idhist+11, NJETC);
	
	// check partons
	Int32_t IQUAC = 0, ICJET = 0;
	for (int i=6; i<parts.size(); ++i) {
		const Particle& part = parts[i];
		
		if (isHardProcess(parts, i) 
		&& part.type == PT_CJET) {
			PT = part.pT();
			ETA = part.getEta(); 
			PHI = part.getPhi();

			if (abs(ETA) < ETCMAX && PT > ETJET) {
				IQUAC++;
				DR = 18.0;

				for (int j=0; j<orecord.Jets.size(); ++j) {
					if (orecord.Jets[j].type == PT_CJET) {  // skoro iterujemy sie po CJetach to po co sprawdzac to jeszcze raz
						DDR = sqrt(
							pow(ETA - orecord.Jets[j].eta_rec, 2) +
							pow(PHI - orecord.Jets[j].phi_rec, 2)
						);
													
						if (abs(PHI - orecord.Jets[j].phi_rec) > PI)
							DDR = sqrt(
								pow(ETA - orecord.Jets[j].eta_rec, 2) +
								pow(abs(PHI - orecord.Jets[j].phi_rec) - 2*PI, 2)
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

	// fill histogram
	histoManager
		->insert(idhist+21, IQUAC );

	for (int i=0; i<orecord.Jets.size(); ++i) {
		PTREC = 0;
		DETRMIN = RCONE;

		if (orecord.Jets[i].type == PT_CJET) {
			for (int j=6; j<parts.size(); ++j) {
				const Particle& part = parts[j];
				
				if (isHardProcess(parts, j) 
				|| part.type != PT_CJET) 
					continue;
				
				PT = part.pT();
				ETA = part.getEta();
				PHI = part.getPhi();

				DPHIA = abs(orecord.Jets[i].phi - PHI); 
				if (DPHIA > PI) 
					DPHIA -= 2*PI;

				DETR = sqrt(
					pow(ETA - orecord.Jets[i].eta , 2) +
					pow(DPHIA, 2)
				);
				
				if (DETR > DETRMIN) 
					continue;

				PTREC = PT;
				DETRMIN = DETR;
			}
		}

		if (PTREC) {
			histoManager
				->insert(idhist + 23, DETRMIN, 1.0);
			histoManager
				->insert(idhist + 24, orecord.Jets[i].pT / PTREC, 1.0);
		}
	}
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
	}
}
