#include "Electron.h"
#include <cstdio>

using namespace AcerDet::analyse;

Electron::Electron( const Configuration& config ) :
	ETCLU	( config.Cluster.MinEt ),
	RCONE	( config.Cluster.ConeR ),

	PTLMIN	( config.Electron.MinMomenta ),
	ETAMAX	( config.Electron.MaxEta ),
	RJE	( config.Electron.MinJetsRlj ),
	RISOLJ	( config.Electron.MinIsolRlj ),
	RDEP	( config.Electron.ConeR ),
	EDMAX	( config.Electron.MaxEnergy ),

	KEYHID	( config.Flag.HistogramId ),
	KEYSME	( config.Flag.Smearing ),

	IEVENT	( 0 ),
	
	histo_isol	("Electron: multiplicity isolated", 0.0, 10.0, 10),
	histo_hard	("Electron: multiplicity hard", 0.0, 10.0, 10),
	histo_sum	("Electron: multipliciyt isol+hard", 0.0, 10.0, 10)
{}

Electron::~Electron() {}

void Electron::printInfo() const {
	// print out title
	printf ("*****************************************\n");
	printf ("*                                       *\n");
	printf ("*     *****************************     *\n");
	printf ("*     ***   analyse::Electron   ***     *\n");
	printf ("*     *****************************     *\n");
	printf ("*                                       *\n");
	printf ("*****************************************\n");

	// print out basic params
	printf ("\n\t... electron isolation ...\n");
	printf ("min. lepton p_T %lf\n", PTLMIN);
	printf ("max. lepton eta %lf\n", ETAMAX);
	printf ("max R_ej for ele-clu %lf\n", RJE);
	printf ("min R_lj for isolation %lf\n", RISOLJ);
	printf ("R for energy deposit %lf\n", RDEP);
	printf ("max E_dep for isolation %lf\n", EDMAX);
	printf ("smearing %s\n", KEYSME ? "on" : "off");
	printf ("\n");

}

void Electron::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	//printf ("Electron: analyse record\n");
	/*
      SUBROUTINE ACDELE(MODE,LPAR,YPAR)

      INTEGER I, II, III, NSTART, NSTOP, IELE, IELEISO
      REAL  DR, DDR, EDEP
      REAL ANGLE
      REAL ETA, PHI, PT, PTCRU
      REAL RJE, RISOLJ, RDEP, EDMAX, PTLMIN, ETAMAX, ETCLU
      REAL RCONE
      REAL ENER, JPT, JPHI, JETA
      REAL PXELE, PYELE, PZELE, EEELE
      REAL RESELE, SIGPH, ENE
      INTEGER KDUM,IDU,IMU, IMAX
      REAL PDUM, ETMAX
      INTEGER LCLU
      LOGICAL ISOL
      INTEGER NBINA, IDENT
      REAL TMAXA,TMINA

// Analiza przypadku
      ELSEIF(MODE.EQ.0.OR.MODE.EQ.2) THEN
C     ===================================
      IEVENT = IEVENT+1


      NELE=0     // liczba znalezionych elektronow w tym obiegu
      DO I=1,100
        KELE(I,2)=0  // macierz K - na flagi    P - na zmienne momenty pedy energie (real)
      ENDDO

// wyznacz koniec historii i poczatek fianl state particle
// analiza dwoch kawalkow eventu osobno
	NSTOP = 0  
	NSTART = 1
	DO I from 1 to N {
		IF(K(I,1) != 21) { 				// K(1) == 21 -> status historia
			NSTOP = I - 1
			NSTART = I
			break
		}
	}

c....look for isolated electrons, sort clusters common
// w czesci drugiej
// K(1) status czastki
	DO I from NSTART to N {
		IF(K(I,1) < 0 || K(I,1) > 10) 
			continue   				// tylko dla czastek ze statusem 1 do 10

c -- analyse electrons
		IF(ABS(K(I,2)) == 11) {    			// K(2) typ czastki +- 11
			ISOL  = TRUE  		        	// izolowany
			LCLU  = 0
			PT = SQRT(P(I,1)**2+P(I,2)**2)   	// zmienne kinematyczne
        		ETA = SIGN(LOG((SQRT(PT**2+P(I,3)**2)+ABS(P(I,3)))/PT),P(I,3))
        		PHI = ANGLE(P(I,1),P(I,2))
        		PTCRU = PT  				// wartosc pedu poprzecznego przed rozmazywaniem
c.....apply smearing
			IF(KEYSME) {          			// jesli ustawiona flaga  rozmaz wg gaussa
				ENE = P(I,4)
				IF(ENE <= 0.0) 
					contninue

				SIGPH = RESELE(ENE,PT,ETA,PHI)
				PXELE = P(I,1) * (1.+SIGPH)
				PYELE = P(I,2) * (1.+SIGPH)
				PZELE = P(I,3) * (1.+SIGPH)
				EEELE = P(I,4) * (1.+SIGPH)
				PT  = SQRT(PXELE**2+PYELE**2)    	// jescze raz wylicz
			}
        
			IF(PT < PTLMIN)              			// stala z konfiga
				contninue
			IF(ABS(ETA) > ETAMAX)			        // stala z konf
				contninue

c === mark eletron-cluster              			        // klastry powiny byc wczesniej zrekonstruowane (powinien byc mapping kluster -elektron)
			DR = 100.0					// znajdz kluster odpowiadajacy elektronowi jako rozwartosc katowa  PCLU - zrekonstruowane klastry
			DO II from 1 to NCLU {				// klastry odpowiadaja rejonom kalorymetry
				DDR = SQRT((ETA-PCLU(II,3))**2+(PHI-PCLU(II,4))**2)
				IF(ABS(PHI-PCLU(II,4)) > PI)
					DDR = SQRT(  (ETA-PCLU(II,3))**2 + (ABS(PHI-PCLU(II,4))-2*PI)**2 )
				IF(DDR < DR) LCLU = II
				DR = MIN(DDR,DR)			// minimalna odleglosc
			}

			IF(DR > RJE) 
				LCLU = 0 				// jesli nie ma klustra w granicach konfiguracji to oznacz jako poza ( nie dotarl do kalorymetru )

c === check if electron  isolated from clusters
			DO II from 1 to NCLU {				// szukaj drugiego najblizszego i sprawdz czy sie miesci w granicach
				IF(II == LCLU) {
					DDR = 100.0
				} ELSE {
					DDR = SQRT((ETA-PCLU(II,3))**2+(PHI-PCLU(II,4))**2)
					IF(ABS(PHI-PCLU(II,4)) > PI)
						DDR = SQRT(  (ETA-PCLU(II,3))**2 + (ABS(PHI-PCLU(II,4))-2*PI)**2 )
				}

				IF(DDR < RISOLJ) 
					ISOL = FALSE
			}

c ===   check on energy deposition of   cells EDEP in cone RDEP		// kluster to grupka celli  cell to fragment kalorymetru
			EDEP = 0.0					// szukaj celli nie uzytych do klastrow
			DO II from 1 to NCELL {
				DDR = SQRT((ETA-PCELL(II,1))**2+(PHI-PCELL(II,2))**2)
				IF(ABS(PHI-PCELL(II,2)) > PI)
					DDR = SQRT(  (ETA-PCELL(II,1))**2 + (ABS(PHI-PCELL(II,2))-2*PI)**2 )
				IF(DDR < RDEP) 
					EDEP += PCELL(II,5)
			}

			IF(EDEP-PTCRU > EDMAX) 
				ISOL = FALSE

c ===  fill /ISOELE/ with isolated electron 
c ===  remove ele-cluster from /CLUSTER/
			IF(ISOL) {						// jest izolowany i od clustrow i od cell wiec laduje na kolekcji
				NELE = NELE + 1
				IF(LCLU != 0) { 			// z listy klustrow usuwamy ten ktory byl z nim stowarzyszony
					DO II from LCLU to NCLU-1 {			// i zapamietaj jego
						DO III from 1 to 5 {
							KCLU(II,III) = KCLU(II+1,III)       // usun i przesun tablice nadpisujac 
                   					PCLU(II,III) = PCLU(II+1,III)
                 				}
					}
					NCLU = NCLU-1
				}

				KELE(NELE,1) = NELE
				KELE(NELE,2) = K(I,2)			// typ elktronu
				KELE(NELE,3) = I
				KELE(NELE,4) = K(K(I,3),2)		// typ matki (w drzewku)
				KELE(NELE,5) = 1

				PELE(NELE,1) = ETA
				PELE(NELE,2) = PHI
           			PELE(NELE,3) = ETA
           			PELE(NELE,4) = PHI
           			PELE(NELE,5) = PT
			}
		}
	}

// dane do histogramu
	CALL HF1(IDENT+11, REAL(NELE), 1.0)

c.....arrange electrons in falling E_T sequence			// sort by E_T
	DO IDU from 1 to NELE {					// naiwne sortowanie przez wybor maxa, wykorzystuje poczatkowy fragment dum jako tmp tablice
		ETMAX = 0
		DO IMU from 1 to NELE {
			IF(KELE(IMU,5) == 0) continue
			IF(PELE(IMU,5) < ETMAX) continue
			IMAX = IMU
			ETMAX = PELE(IMU,5)
		}

		KELE(IMAX,5) = 0
		DO II from 1 to 5 {
			KDUM(IDU,II) = KELE(IMAX,II)
			PDUM(IDU,II) = PELE(IMAX,II)
		}
	}

	DO I from 1 to NELE {
		DO II from 1 to 5 {				// skopiuj dum -> ele
			KELE(I,II) = KDUM(I,II)
			PELE(I,II) = PDUM(I,II)
		}
		KELE(I,1) = I					// numer elektrona
		KELE(I,5) = 1					// stan(1) = przetworzony?
	}


// dla czastek z histori
	IELE    = 0
	IELEISO = 0

	DO I from 1 to NSTOP {
		IF(ABS(K(I,2)) == 11) {
			PT = SQRT(P(I,1)**2+P(I,2)**2) 
			ETA = SIGN(LOG((SQRT(PT**2+P(I,3)**2)+ABS(P(I,3)))/PT),P(I,3)) 
			PHI = ANGLE(P(I,1),P(I,2))
			ENER = 0.0
			ISOL = TRUE

			DO II from 1 to NSTOP {
				IF(ABS(K(II,2)) <= 21 && I != II && ABS(K(II,2)) != 12 && ABS(K(II,2)) != 14 && ABS(K(II,2)) != 16) {
					// przelicz dane dla II czastki - standard per particle
					JPT = SQRT(P(II,1)**2+P(II,2)**2) 
					JETA = SIGN(LOG((SQRT(JPT**2+P(II,3)**2)+ABS(P(II,3)))/JPT),P(II,3)) 
					JPHI = ANGLE(P(II,1),P(II,2))
					DDR = SQRT((ETA-JETA)**2+(JPHI-PHI)**2)

					IF(ABS(JPHI-PHI) > PI)
						DDR = SQRT((ETA-JETA)**2 + (ABS(JPHI-PHI)-2*PI)**2 )

					IF(DDR < RISOLJ && JPT > ETCLU) 
						ISOL = FALSE

					IF(DDR < RDEP && JPT < ETCLU) 
						ENER += JPT
				}
			}

			IF(ENER > EDMAX) 
				ISOL = FALSE

			IF(ABS(ETA) < ETAMAX && PT > PTLMIN)
				IELE = IELE + 1

			IF(ABS(ETA) < ETAMAX && PT > PTLMIN && ISOL)
				IELEISO = IELEISO + 1
		}
	}

// ile takich ele bylo w tej czesci i ile bylo izolowanych (check na tym co robilimsy wczesniej) 
         CALL HF1(IDENT+21, REAL(IELE), 1.0)
         CALL HF1(IDENT+31, REAL(IELEISO), 1.0)

C
      ELSEIF(MODE.EQ.1) THEN
C     =======================

// druk histogramow

      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '        OUTPUT FROM AcerDET      '
      WRITE(NOUT,BXTXT) '              ACDELE             '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXCLO)

      CALL HPRINT(IDENT+11)
      CALL HPRINT(IDENT+21)
      CALL HPRINT(IDENT+31)

      ENDIF
C     =====

      END
*/

}
