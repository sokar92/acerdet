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
	//printf ("CJet: analyse record\n");
	/*
      SUBROUTINE ACDCJE(MODE,LPAR,YPAR)
      
      INTEGER KEYHID, KEYBCL, IEVENT 
      INTEGER I, II, NSTART,NSTOP
      REAL  DR, DDR
      REAL ANGLE
      REAL ETA, PHI, PT
      REAL RJC, PTCMIN, ETCMAX, ETJET
      INTEGER NJETC
      INTEGER JETC
      LOGICAL CJET
      INTEGER IQUAC, ICJET
      REAL RCONE, PTREC,DETRMIN, DETR, DPHIA
      INTEGER IJET
      INTEGER NBINA
      INTEGER IDENT
      REAl TMAXA,TMINA

      ELSEIF(MODE == 0 && KEYBCL == true) THEN
C     ========================================
	IEVENT = IEVENT + 1

	NSTOP = 0
	NSTART = 1
	DO I = from 1 to N {
		IF(K(I,1) != 21) {
			NSTOP = I-1
			NSTART = I
			break
		}
	}

	NJETC = 0

c....look for c-jets
	DO I = from NSTART to N {							// przegladaj nie historyczne
		IF(ABS(K(I,2)) == 4 && K(I,1) != 21) {					// if (p->type() == C_JET && p->status() != PS_HISTORY) {

c....if there is a c-quark found before hadronization
c....if there are still jets
			IF(NJET > 0) {
				CJET = TRUE
				JETC = 0

c....and this c-quark is the last one in the FSR cascade
				IF (K(I,4) != 0) {							// istnieje corka
					DO II = form K(I,4) to K(I,5) {					// iteruj po corkach
						IF (ABS(K(II,2)) == 4) 
							CJET = FALSE					// jakas corka jets cjetem
					}
				}

				IF(!CJET) continue

				PT = SQRT(P(I,1)*2+P(I,2)*2) 						// pt = srqt( px^2 + py^2 )
				IF(PT < PTCMIN) CJET = FALSE						// pt mniejsze niz minimum z configa
				IF(!CJET) continue

				ETA = SIGN(LOG((SQRT(PT*2+P(I,3)*2)+ABS(P(I,3)))/PT),P(I,3)) 
				IF(ABS(ETA) > ETCMAX) CJET = FALSE					// kat poza zakresem stozka z configa
				IF(!CJET) continue

				PHI = ANGLE(P(I,1),P(I,2))						// kat pomiedzy px i py

c === mark c-jet
				DR = 100.0
				DO II = from 1 to NJET {
					IF(ABS(KJET(II,2)) != 5) {					// ma sens ze wzgledu na kolejnosc B->C  if (p->type() != B_JET) {
						DDR = SQRT((ETA-PJET(II,3))*2+(PHI-PJET(II,4))*2)
						IF(ABS(PHI-PJET(II,4)) > PI)
							DDR = SQRT(  (ETA-PJET(II,3))*2+(ABS(PHI-PJET(II,4))-2*PI)*2 )
						IF(DDR < DR) 
							JETC = II
						DR = MIN(DDR,DR)
					}
				}

				IF(DR > RJC) {
					CJET = FALSE
					JETC = 0
				}
				IF(!CJET) continue

c ===  labell  c-jet
				KJET(JETC,2) = 4
				KJET(JETC,5) = I
				NJETC = NJETC + 1							// kolejny cjet znaleziony
			}
		}
	}

	CALL HF1(IDENT+11,REAL(NJETC),1.0)					// zapis do histogramu

c === check partons
	IQUAC = 0
	ICJET = 0
	DO I = from 7 to NSTOP {						// CZEMU OD 7 ?
		IF(ABS(K(I,2)) == 4) {						// if (p->type() == C_JET) {
			PT = SQRT(P(I,1)*2+P(I,2)*2)						// pt = sqrt( px^2 + py^2 ) 
			ETA = SIGN(LOG((SQRT(PT*2+P(I,3)*2)+ABS(P(I,3)))/PT),P(I,3)) 
			PHI = ANGLE(P(I,1),P(I,2))						// kat pomiedzy px i py

			IF(ABS(ETA) < ETCMAX && PT > ETJET) {					// spelnia widelki z configa
				IQUAC = IQUAC + 1
				DR = 18.0

				DO II = from 1 to NJETC {
					IF(ABS(KJET(II,2)) == 4) {
						DDR = SQRT((ETA-PJET(II,3))*2+(PHI-PJET(II,4))*2)
						IF(ABS(PHI-PJET(II,4)) > PI)
							DDR = SQRT( (ETA-PJET(II,3))*2+(ABS(PHI-PJET(II,4))-2*PI)*2 )
						IF(DDR < DR) 
							ICJET = II
						IF(DDR < DR) 
							DR = DDR
					}
				}
			}
		}
	}

	CALL HF1(IDENT+21,REAL(IQUAC),1.0)					// zapis do histogramu

	DO IJET = from 1 to NJET {
		PTREC = 0
		DETRMIN = RCONE

		IF(ABS(KJET(IJET,2)) == 4) {						// if (p->type() == C_JET)
			DO I = from 7 to N { 						// for (int i=7;i<n;++i)
				IF(K(I,1) != 21 || ABS(K(I,2)) != 4) continue		// if (p->status() != PS_HISTORY || p->type() != C_JET)
				PT = SQRT(P(I,1)**2+P(I,2)**2)						// pt = sqrt( px^2 + py^2 )       ** = ^ ?					
				ETA = SIGN(LOG((SQRT(PT*2+P(I,3)*2)+ABS(P(I,3)))/PT),P(I,3)) 
				PHI = ANGLE(P(I,1),P(I,2))						// kat pomiedzy px i py

				DPHIA = ABS(PJET(IJET,4)-PHI) 
				IF(DPHIA > PI) 
					DPHIA = DPHIA-2*PI

				DETR = SQRT((ETA-PJET(IJET,3))**2+DPHIA**2)
				IF(DETR > DETRMIN) continue

				PTREC = PT
				DETRMIN = DETR
			}
		}

		IF(PTREC != 0) {
			CALL HF1(IDENT+23,DETRMIN,1.0)
			CALL HF1(IDENT+24,PJET(IJET,5)/PTREC,1.0)
		}
	}

C
      ELSEIF(MODE == 1 && KEYBCL != 0) THEN
C     ========================================

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
