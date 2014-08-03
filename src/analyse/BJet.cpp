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
	/*
	const vector<Particle>& parts = irecord.particles();
	int n = parts.size();
	
	int nstart = 0;
	while(nstart < n && parts[nstart].state != PS_HISTORY) 
		nstart++;
	int nstop = nstart-1;*/
	
	/*
      SUBROUTINE ACDBJE(MODE,LPAR,YPAR)
      
      INTEGER I, II, NSTART,NSTOP
      REAL  DR, DDR
      REAL ANGLE
      REAL ETA, PHI, PT
      REAL PTBMIN, ETBMAX, RJB, ETJET
      INTEGER NJETB
      INTEGER JETB
      LOGICAL BJET
      INTEGER IQUAB, IBJET
      REAL RCONE, PTREC,DETRMIN, DETR, DPHIA
      INTEGER IJET
      INTEGER NBINA
      INTEGER IDENT
      REAl TMAXA,TMINA

processing ->
 
->      NSTOP=0
->      NSTART=1
->      DO I=1, N
->       IF(K(I,1) != 21) THEN
->           NSTOP = I-1
->           NSTART= I
->           GOTO 500
->       ENDIF
->      ENDDO
-> 500  CONTINUE

        NJETB = 0								// na poczatku jest 0 b-jetow

c....look for b-jets
	DO I = from NSTART to N {						// przegladaj nie historyczne elementy ' for (int i=start; i<n;++i) { '
		IF(ABS(K(I,2)) == 5 && K(I,1) != 21) {				// K(2) = 5 | -5 -> id czastki	' if (p->type() == B_JETS && p->status() != PS_HISTORY) { '

c....if there is a b-quark found before hadronization
c....if there are still jets
			IF(NJET > 0) {
				BJET = TRUE
				JETB = 0

c....and this b-quark is the last one in the FSR cascade
				IF (K(I,4) != 0) {					// istnieje corka w drzewku
					DO II = from K(I,4) to K(I,5) {			// przegladaj corki
						IF (ABS(K(II,2)) == 5) 
							BJET = FALSE			// jesli typ ktorejkolwiek to B_JET to glowny nie jest B_JETEM
					}
				}

				IF(!BJET) continue					// jesli nie jest B_JETEM to koniec

				PT = SQRT(P(I,1) * 2 + P(I,2) * 2) 
				IF(PT < PTBMIN) BJET = FALSE  					// ponizej normy (z config)
				IF(!BJET) continue						// nie jest

				ETA = SIGN(LOG((SQRT(PT*2+P(I,3)*2)+ABS(P(I,3)))/PT),P(I,3)) 	// wyznacz kat	
				IF(ABS(ETA) > ETBMAX) BJET = FALSE					// kat poza zakresem (z config)
				IF(!BJET) continue							// nie jest
			
				PHI = ANGLE(P(I,1),P(I,2))						// kat pomiedzy px i py

c === mark b-jet
				DR = 100.0
				DO II = from 1 to NJET {
					DDR = SQRT((ETA-PJET(II,3))*2+(PHI-PJET(II,4))*2)
					IF( ABS(PHI-PJET(II,4)) > PI )
						DDR = SQRT( (ETA-PJET(II,3))*2 + (ABS(PHI-PJET(II,4))-2*PI)*2 )		// jest na to juz procedurka
					IF(DDR < DR) 
						JETB = II
					DR = MIN(DDR,DR)
				}

				IF(DR > RJB) {
					BJET = FALSE
					JETB = 0
				}

				IF(!BJET) continue							// nie jest b-jetem

c ===  labell  b-jet
				KJET(JETB,2) = 5
				KJET(JETB,5) = I
				NJETB = NJETB + 1							// jest b-jetem -> zwieksz licznik wykrytych
			}
		}
	}

	CALL HF1(IDENT+11, REAL(NJETB), 1.0)							// zapis do histogramow

c === check partons
	IQUAB = 0
	IBJET = 0
	DO I = from 7 to NSTOP {								// CZEMU OD 7? -> barcode = numery czastek w zdarzeniu
		IF(ABS(K(I,2)) == 5) {								// jest to B_JET -> przetworz
			PT = SQRT(P(I,1)*2+P(I,2)*2)						// pt = sqrt( px^2 + py^2 ) 
			ETA = SIGN(LOG((SQRT(PT**2+P(I,3)**2)+ABS(P(I,3)))/PT),P(I,3)) 
			PHI = ANGLE(P(I,1),P(I,2))						// kat pomiedzy px i py

			IF(ABS(ETA) < ETBMAX && PT > ETJET) {					// miesci sie w granicach z configa
				IQUAB = IQUAB + 1
				DR = 18.0

				DO II = from 1 to NJET {
					IF(ABS(KJET(II,2)) == 5) {				// jest b_jetem
						DDR = SQRT((ETA-PJET(II,3))*2+(PHI-PJET(II,4))*2)		// jest na to procedura
						IF( ABS(PHI-PJET(II,4)) > PI )
							DDR = SQRT(  (ETA-PJET(II,3))*2 +(ABS(PHI-PJET(II,4))-2*PI)*2 )
						IF(DDR < DR) 
							IBJET = II
						IF(DDR < DR) 
							DR = DDR
					}
				}
			}
		}
	}

	CALL HF1(IDENT+21,REAL(IQUAB),1.0)					// zapis do histo

	DO IJET = from 1 to NJET {
		PTREC = 0
		DETRMIN = RCONE

		IF(ABS(KJET(IJET,2)) == 5) {						// jest b_jetem
			DO I = from 7 to N {						// CZEMU OD 7 
				IF(K(I,1) != 21 || ABS(K(I,2)) != 5) continue		// if (p->status() != 21 || p->type() != B_JET) {

				PT = SQRT(P(I,1)*2+P(I,2)*2)						// pt = sqrt( px^2 + py^2 )
				ETA = SIGN(LOG((SQRT(PT*2+P(I,3)*2)+ABS(P(I,3)))/PT),P(I,3)) 
				PHI = ANGLE(P(I,1),P(I,2))						// kat pomiedzy px i py
				DPHIA = ABS(PJET(IJET,4)-PHI) 

				IF(DPHIA > PI) 
					DPHIA = DPHIA-2*PI					// kat w przedziale [-pi,pi]

				DETR = SQRT((ETA-PJET(IJET,3))*2+DPHIA*2)
				IF(DETR > DETRMIN) continue 

				PTREC = PT
				DETRMIN = DETR
			}
		}

		IF( PTREC != 0 ) {
			CALL HF1(IDENT+23,DETRMIN,1.0)
			CALL HF1(IDENT+24,PJET(IJET,5)/PTREC,1.0)
		}
	}

end ->
      ELSEIF(MODE.EQ.1.AND.KEYBCL.NE.0) THEN

      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '        OUTPUT FROM AcerDET      '
      WRITE(NOUT,BXTXT) '              ACDBJE             '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXCLO)


      CALL HPRINT(IDENT+11)
      CALL HPRINT(IDENT+21)

      ENDIF
*/

}

void BJet::printResults() const {
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
