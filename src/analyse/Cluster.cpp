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
	Real64_t PTLRAT = 1.0 / pow(sinh(ETACLU), 2.0);
	
	// Find initiator cell: the one with highest pT of not yet used ones.
	Real64_t ETMAX = 0.0;
	Int32_t ICMAX = -1;
	Real64_t ETA, PHI, THETA, PT;
	
	for (int i=0; i<orecord.vCell.size(); ++i) {
		if (orecord.vCell[i].K[4] != 2) continue;
		if (orecord.vCell[i].P[4] > ETMAX) {
			ICMAX = i;
			ETA = orecord.vCell[i].P[0];
			PHI = orecord.vCell[i].P[1];
			ETMAX = orecord.vCell[i].P[4];
		}
	}
	
	// maximum found Et is not enough
	if (ETMAX < ETINI) {
	}
	
	/*
       IEVENT = IEVENT + 1

->      PTLRAT = 1. / SINH(ETACLU)**2 

C...Find initiator cell: the one with highest pT of not yet used ones. 
      
->	NDUM = 0 
->  	150 ETMAX = 0. 								// jakis DO ... WHILE ... 
->	DO IC = from 1 to NCELL {							// przeszukuj wszystkie celki 
->		IF(KCELL(IC,5) != 2) continue;
->		IF(PCELL(IC,5) <= ETMAX) continue; 
->		ICMAX = IC 								// ostatni numer cellki nie spelniajacy warunkow
->		ETA = PCELL(IC,1)							// dane wez z tej celki 
->		PHI = PCELL(IC,2) 
->		ETMAX = PCELL(IC,5) 
->	}

->      IF(ETMAX < ETINI) GOTO 220					// etmax jest zbyt male (ponizej init limitu) 

	KCELL(ICMAX,5) = 1 
	NDUM = NDUM + 1 
	KDUM(NDUM,3) = 1
	KDUM(NDUM,4) = KCELL(ICMAX,4) 
	KDUM(NDUM,5) = 1 
	
	PDUM(NDUM,1) = ETA 
	PDUM(NDUM,2) = PHI 
	PDUM(NDUM,3) = 0 
	PDUM(NDUM,4) = 0 
	PDUM(NDUM,5) = 0 
 
C...Sum up unused cells within required distance of initiator. 
	DO IC = from 1 to NCELL {						// przegladaj po celkach i wykonuj pewne aktualizacje 
		IF(KCELL(IC,5) == 0) continue 
		DPHIA = ABS(PCELL(IC,2)-PHI) 
		PHIC = PCELL(IC,2) 

		IF(DPHIA > PI) 
			PHIC = PHIC + SIGN(TWOPI,PHI)

		IF(ABS(ETA) < CALOTH && (PCELL(IC,1)-ETA)**2+(PHIC-PHI)**2.GT.RCONE**2) continue
		IF(ABS(ETA) > CALOTH && (PCELL(IC,1)-ETA)**2+(PHIC-PHI)**2.GT.RCONE**2) continue 

		KCELL(IC,5) = -KCELL(IC,5) 
		KDUM(NDUM,3) = KDUM(NDUM,3) + 1 
		KDUM(NDUM,4) = KDUM(NDUM,4) + KCELL(IC,4)
		PDUM(NDUM,3) = PDUM(NDUM,3) + PCELL(IC,5) * PCELL(IC,1) 
		PDUM(NDUM,4) = PDUM(NDUM,4) + PCELL(IC,5) * PHIC 
		PDUM(NDUM,5) = PDUM(NDUM,5) + PCELL(IC,5) 
	}

 
C...Reject cluster below minimum ET, else accept. 
	IF(PDUM(NDUM,5) < ETCLU) 
	THEN 
		NDUM = NDUM - 1 
		DO IC = from 1 to NCELL	{						// lec po celkach  
			IF(KCELL(IC,5) < 0) 
				KCELL(IC,5) = -KCELL(IC,5)				// a = |a| 
		}
	ELSE
		PDUM(NDUM,3) = PDUM(NDUM,3) / PDUM(NDUM,5) 
		PDUM(NDUM,4) = PDUM(NDUM,4) / PDUM(NDUM,5) 

		IF(ABS(PDUM(NDUM,4)) > PI) 
			PDUM(NDUM,4) = PDUM(NDUM,4) - SIGN(TWOPI,PDUM(NDUM,4)) 

		DO IC = from 1 to NCELL {						// petla po celkach  
C>>>>>> mark used cells in /CELLS/
			IF(KCELL(IC,5) < 0) 
				KCELL(IC,5) = 0						
		}
	ENDIF 
      GOTO 150 

 220  CONTINUE
c....Arrange clusters in falling ET sequence and store clusters in /CLUSTER/
	NCLU = 0
	DO ICLU = from 1 to NDUM {
		ETMAX = 0
		DO IDUM = from 1 to NDUM {						// szukaj maxa
			IF(KDUM(IDUM,5) == 0) continue
			IF(PDUM(IDUM,5) < ETMAX) continue
			IMAX = IDUM
			ETMAX = PDUM(IDUM,5)
		}
		KDUM(IMAX,5) = 0 //sorted state
		NCLU = NCLU + 1								// zapisz nowy cluster na outputrecord
		DO II = from 1 to 5 {
			KCLU(NCLU,II) = KDUM(IMAX,II)
			PCLU(NCLU,II) = PDUM(IMAX,II)
		}
		KCLU(NCLU,5) = 1
	}

c.....histogram NCLU
      CALL HF1(IDENT+1,REAL(NCLU),1.0)						// zapis do histogramu

c....reconstruct baricenter of particles
	DO ICLU = from 1 to NCLU {							// petla po wszytskich clustrach
		ETAREC = 0
		PTREC = 0
		PHIREC = 0

		DO I = from 1 to N {
			IF(K(I,1) <= 0 || K(I,1) > 10) continue 
			IF(P(I,1)**2+P(I,2)**2 <= PTLRAT*P(I,3)**2) contninue
 
			KC = KFCOMP(K(I,2)) 
			IF(KC == 0 || KC == 12 || KC == 14 || KC == 16 || KC == 13 || KC == KFINVS) continue

			DETPHI = 0.0
			IF(KEYFLD != false && KUCHGE(K(I,2)) != 0) 
			THEN									// flaga z configa
				PT = SQRT(P(I,1)**2+P(I,2)**2)					// pt = sqrt( px^2 + py^2 )
				IF(PT < PTMIN) contninue					// continue

				ETA = SIGN(LOG((SQRT(PT**2+P(I,3)**2)+ABS(P(I,3)))/PT),P(I,3)) 
				PHI = ANGLE(P(I,1),P(I,2))					// kat pomiedzy px i py
				CHRG = KUCHGE(K(I,2))/3.
				
				IF(ABS(P(I,3)/P(I,4)) < 1) 
				THEN
					THETA = ACOS(P(I,3) / SQRT(PT**2+P(I,3)**2))		// kat z iloczynu skalarnego (pz,p)
				ELSE
					IF(P(I,3) > 0) THETA = 0				// normalizacja do [0,pi]
					IF(P(I,3) < 0) THETA = PI
				ENDIF
				
				DETPHI = CHRG * FLDPHI(PT,ETA,THETA)
			ENDIF 

			PT = SQRT(P(I,1)**2+P(I,2)**2)						// pt = sqrt( px^2 + py^2 )
			ETA = SIGN(LOG((SQRT(PT**2+P(I,3)**2)+ABS(P(I,3)))/PT),P(I,3)) 
			PHI = ANGLE(P(I,1),P(I,2))						// kat pomiedzy px i py
			PHI = PHI + DETPHI
			DPHIA = ABS(PCLU(ICLU,2)-PHI) 

			IF(DPHIA > PI) 
				DPHIA = DPHIA-2*PI					// normalizacja do [-pi,pi]

			IF(ABS(PCLU(ICLU,1)) < CALOTH && (PCLU(ICLU,1)-ETA)**2+(DPHIA)**2 > RCONE**2) contninue
			IF(ABS(PCLU(ICLU,1)) > CALOTH && (PCLU(ICLU,1)-ETA)**2+(DPHIA)**2 > RCONE**2) continue

			PTREC = PTREC + PT
			ETAREC = ETAREC + ETA * PT
			PHIREC = PHIREC + PHI * PT
		}

		ETAREC /= PTREC
		PHIREC /= PTREC
		DETR = SQRT((ETAREC-PCLU(ICLU,3))**2+(PHIREC-PCLU(ICLU,4))**2)
      
		CALL HF1(IDENT+11, ETAREC-PCLU(ICLU,3), 1.0)				// zapis do histogramow
		CALL HF1(IDENT+12, PHIREC-PCLU(ICLU,4), 1.0)
		CALL HF1(IDENT+13, DETR, 1.0)
		CALL HF1(IDENT+14, PCLU(ICLU,5)/PTREC, 1.0)
	}

	DO ICLU from 1 to NCLU { 						// lec po wszytskich klastrach
		PTREC = 0
		DETRMIN = RCONE

		DO I from 7 to N { 
			IF(K(I,1) != 21 || ABS(K(I,2)) > 10) contninue

			PT = SQRT(P(I,1)**2+P(I,2)**2)
			ETA = SIGN(LOG((SQRT(PT**2+P(I,3)**2)+ABS(P(I,3)))/PT),P(I,3)) 
			PHI = ANGLE(P(I,1),P(I,2))
			DPHIA = ABS(PCLU(ICLU,4)-PHI) 
			IF(DPHIA > PI) 
				DPHIA = DPHIA-2*PI
			DETR = SQRT((ETA-PCLU(ICLU,3))**2+DPHIA**2)
			IF(DETR > DETRMIN) continue 
			PTREC = PT
			DETRMIN = DETR
		}

		IF(PTREC) 
		THEN
			CALL HF1(IDENT+23,DETRMIN,1.0)
			CALL HF1(IDENT+24,PCLU(ICLU,5)/PTREC,1.0)
		ENDIF
	}

	 */
}
