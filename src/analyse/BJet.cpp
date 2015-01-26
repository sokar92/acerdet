#include "BJet.h"
#include <cstdio>
using namespace AcerDet::analyse;

#include "../core/Functions.h"
using namespace AcerDet::core;

BJet::BJet( const Configuration& config, IHistogramManager *histoMng ) :
	ETJET	( config.Jet.MinEnergy ),
	RCONE	( config.Cluster.ConeR ),

	PTBMIN	( config.BJet.MinMomenta ),
	ETBMAX	( config.BJet.MaxEta ),
	RJB	( config.BJet.MaxRbj ),

	KEYHID	( config.Flag.HistogramId ),
	KEYBCL	( config.Flag.BCJetsLabeling ),

	IEVENT	( 0 ),

	histoManager( histoMng ),
	histoRegistered( false )
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
	
	Int32_t idhist = 700 + KEYHID;
	if (!histoRegistered) {
		histoRegistered = true;
		histoManager
			->registerHistogram(idhist+11, "BJet: b-jets multiplicity", 10, 0.0, 10.0);
		histoManager
			->registerHistogram(idhist+21, "BJet: b-quarks HARD multiplicity", 10, 0.0, 10.0);
		histoManager
			->registerHistogram(idhist+23, "BJet: delta r bjet-cquark", 50, 0.0,  5.0);
		histoManager
			->registerHistogram(idhist+24, "BJet: pTbjet/pTbquark", 50, 0.0,  2.0);
	}
	
	// do not use this algorithm
	if (!KEYBCL)
		return;

	// new event to compute
	IEVENT++;

	// reference to particles container
	const vector<Particle>& parts = irecord.particles();
	
	// look for b-jets
	Int32_t NJETB = 0;
	for (int i=0; i<parts.size(); ++i) {
		const Particle& part = parts[i];
if (part.type == PT_BJET) printf ("Status: %d -> st_id %d\n", part.status, part.statusID);		
		if (part.status == PS_CASCADE_QUARK //abs(part.statusID) >= 30//part.status == PS_DECAYED
		&& part.type == PT_BJET) {
printf ("\nFound B-quark\n");
			// if there is a b-quark found before hadronization
			// if there are still jets
			if (!orecord.Jets.empty()) {
printf ("Not empty record\n");				
				Bool_t BJET = true;
				
				// and this b-quark is the last one in the FSR cascade
				if (part.hasDaughter()) {
					for (int j=part.daughters.first; j<=part.daughters.second; ++j) {
						if (parts[j-1].type == PT_BJET)
							BJET = false;
					}
				}

				if (!BJET) 
					continue;
 printf ("Is b-jet\n");
				if (part.pT() < PTBMIN) 
					continue;

				if (abs(part.getEta()) > ETBMAX)
					continue;
printf ("Matches conditions\n");				
				// mark b-jet
				Real64_t DR = 100.0;
				Int32_t JETB = -1;
				
				for (int j=0; j<orecord.Jets.size(); ++j) {
					
					Real64_t DDR = sqrt(
						pow(part.getEta() - orecord.Jets[j].eta_rec, 2) +
						pow(part.getPhi() - orecord.Jets[j].phi_rec, 2)
					);
					
					if (abs(part.getPhi() - orecord.Jets[j].phi_rec) > PI )
						DDR = sqrt(
							pow(part.getEta() - orecord.Jets[j].eta_rec, 2) + 
							pow(abs(part.getPhi() - orecord.Jets[j].phi_rec) - 2*PI, 2)
						);
						
					if (DDR < DR) { 
						JETB = j;
						DR = DDR;
					}
				}

				// label  b-jet
				if (JETB >= 0 && DR <= RJB) {
					// KJET(JETB,2) = 5; // typ = bjet
					// KJET(JETB,5) = I;
					orecord.Jets[JETB].type = B_JET;
					NJETB++;
					printf ("NJETB++\n");
				}
			}
		}
	}

	// histogram store
	histoManager
		->insert(idhist+11, NJETB, 1.0);
	
	// check partons
	Int32_t IQUAB = 0, IBJET = 0;
	for (int i=6; i<parts.size(); ++i) {
		const Particle& part = parts[i];
		
		if (part.status == PS_HP_QUARK //isHardProcess(parts, i)
		&& part.type == PT_BJET) {

			if (abs(part.getEta()) < ETBMAX
			&& part.pT() > ETJET) {
				
				IQUAB++;
				Real64_t DR = 18.0;

				for (int j=0; j<orecord.Jets.size(); ++j) {
					if (orecord.Jets[j].type == B_JET) {
						
						Real64_t DDR = sqrt(
							pow(part.getEta() - orecord.Jets[j].eta_rec, 2) +
							pow(part.getPhi() - orecord.Jets[j].phi_rec, 2)
						);
						
						if (abs(part.getPhi() - orecord.Jets[j].phi_rec) > PI)
							DDR = sqrt(
								pow(part.getEta() - orecord.Jets[j].eta_rec, 2) +
								pow(abs(part.getPhi() - orecord.Jets[j].phi_rec) - 2*PI, 2)
							);
						
						if (DDR < DR) {
							IBJET = j; 
							DR = DDR;
						}
					}
				}
			}
		}
	}

	// fill histogram
	histoManager
		->insert(idhist+21, IQUAB );
	
	for (int i=0; i<orecord.Jets.size(); ++i) {
		Real64_t PTREC = 0.0;
		Real64_t DETRMIN = RCONE;

		if (orecord.Jets[i].type == B_JET) {
			for (int j=6; j<parts.size(); ++j) {
				const Particle& part = parts[j];
				
				if (isHardProcess(parts, j)
				|| part.type != PT_BJET)   // sprawdz czy ma byc || czy &&
					continue;

				Real64_t DPHIA = abs(part.getPhi() - orecord.Jets[i].phi_rec); 
				if (DPHIA > PI) 
					DPHIA -= 2*PI;

				Real64_t DETR = sqrt(
					pow(part.getEta() - orecord.Jets[i].eta_rec, 2) +
					pow(DPHIA, 2)
				);
				
				if (DETR < DETRMIN) {
					PTREC = part.pT();
					DETRMIN = DETR;
				}
			}
		}

		if (PTREC != 0) {
			histoManager
				->insert(idhist + 23, DETRMIN, 1.0);
			histoManager
				->insert(idhist + 24, orecord.Jets[i].pT / PTREC, 1.0);
		}
	}
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
	}
}
