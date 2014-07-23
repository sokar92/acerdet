#include "Photon.h"
#include <cstdio>

using namespace AcerDet::analyse;

Photon::Photon( const Configuration& config ) :
	ETCLU	( config.Cluster.MinEt ),
	RCONE	( config.Cluster.ConeR ),

	PTLMIN	( config.Photon.MinMomenta ),
	ETAMAX	( config.Photon.MaxEta ),
	RJE	( config.Photon.MinJetsRlj ),
	RISOLJ	( config.Photon.MinIsolRlj ),
	RDEP	( config.Photon.ConeR ),
	EDMAX	( config.Photon.MaxEnergy ),

	KEYHID	( config.Flag.HistogramId ),
	KEYSME	( config.Flag.Smearing ),

	IEVENT	( 0 ),
	
	histo_isol	("Photon: multiplicity isolated", 0.0, 10.0, 10),
	histo_hard	("Photon: multiplicity hard", 0.0, 10.0, 10),
	histo_sum	("Photon: multipliciyt isol+hard", 0.0, 10.0, 10)
{}

Photon::~Photon() {}

void Photon::printInfo() const {
	// print out title
	printf ("***************************************\n");
	printf ("*                                     *\n");
	printf ("*     ***************************     *\n");
	printf ("*     ***   analyse::Photon   ***     *\n");
	printf ("*     ***************************     *\n");
	printf ("*                                     *\n");
	printf ("***************************************\n");

	// print out basic params
	printf ("\n\t... photon isolation ...\n");
	printf ("min. photon p_T %lf\n", PTLMIN);
	printf ("max. photon eta %lf\n", ETAMAX);
	printf ("max R_gam-clust %lf\n", RJE);
	printf ("min R_isol %lf\n", RISOLJ);
	printf ("R for energy deposit %lf\n", RDEP);
	printf ("max E_dep for isolation %lf\n", EDMAX);
	printf ("smearing %s\n", KEYSME ? "on" : "off");
	printf ("\n");

}

void Photon::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	//printf ("Photon: analyse record\n");
	/*
      SUBROUTINE ACDPHO(MODE,LPAR,YPAR)

      INTEGER I, II, III, NSTART, NSTOP, IPHO,IPHISO
      REAL  DR, DDR, EDEP
      REAL ANGLE
      REAL ETA, PHI, PT, PTCRU
      REAL RJE, RISOLJ, RDEP, EDMAX, PTLMIN, ETAMAX, ETCLU
      REAL RCONE
      REAL ENER, JPT, JPHI, JETA
      REAl PXPHO, PYPHO, PZPHO, EEPHO,SIGPH
      REAL RESPHO, ENE
      INTEGER KDUM,IDU,IMU, IMAX
      REAL PDUM, ETMAX
      INTEGER LCLU
      LOGICAL ISOL
      INTEGER NBINA, IDENT
      REAL TMAXA,TMINA

      ELSEIF(MODE.EQ.0.OR.MODE.EQ.2) THEN
C     ===================================
	IEVENT = IEVENT + 1

c....look for isolated photons, sort cluster commons
	NPHO = 0
	DO I from 1 to 100 {
        	KPHO(I,2) = 0
	}

	NSTOP = 0
	NSTART = 1
	DO I from 1 to N {
		IF(K(I,1) != 21) {
			NSTOP = I - 1
			NSTART = I
			break
		}
	}

	DO I from NSTART to N {
		IF(K(I,1) < 0 || K(I,1) > 10) continue
		IF(SQRT(P(I,1)**2+P(I,2)**2) == 0) continue

c -- analyse photons
		IF(ABS(K(I,2)) == 22) {
			ISOL  = TRUE
			LCLU  = 0
			PT = SQRT(P(I,1)**2+P(I,2)**2)
			ETA = SIGN(LOG((SQRT(PT**2+P(I,3)**2)+ABS(P(I,3)))/PT),P(I,3))
			PHI = ANGLE(P(I,1),P(I,2))
			PTCRU = PT

c....apply smearing
			IF(KEYSME) {
				ENE = P(I,4)
				IF(ENE <= 0.0) continue
				SIGPH  = RESPHO(ENE,PT,ETA,PHI)
				pxpho  = P(I,1) * (1.+sigph)
				pypho  = P(I,2) * (1.+sigph)
				pzpho  = P(I,3) * (1.+sigph)
				eepho  = P(I,4) * (1.+sigph)
				PT  = sqrt(pxpho**2+pypho**2)
				ENE = eepho
			}

			IF(PT < PTLMIN) 
				contninue 
			IF(ABS(ETA) > ETAMAX) 
				contninue 

c === mark photon-cluster
			DR = 100.0
			DO II from 1 to NCLU {
				DDR = SQRT((ETA-PCLU(II,3))**2+(PHI-PCLU(II,4))**2)
				IF(ABS(PHI-PCLU(II,4)) > PI)
					DDR = SQRT(  (ETA-PCLU(II,3))**2 + (ABS(PHI-PCLU(II,4))-2*PI)**2 )
				IF(DDR < DR) 
					LCLU = II
				DR = MIN(DDR,DR)
			}

			IF(DR > RJE) 
				LCLU = 0
 
c === check if photon  isolated from clusters
			DO II from 1 to NCLU {
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

c ===   check on energy deposition of   cells EDEP in cone RDEP
			EDEP =  0.0
			DO II from 1 to NCELL {
				DDR = SQRT((ETA-PCELL(II,1))**2+(PHI-PCELL(II,2))**2)
				IF(ABS(PHI-PCELL(II,2)) > PI)
					DDR = SQRT(  (ETA-PCELL(II,1))**2 + (ABS(PHI-PCELL(II,2))-2*PI)**2 )
				IF(DDR < RDEP) 
					EDEP += PCELL(II,5)
			}

			IF(EDEP-PT > EDMAX) 
				ISOL = FALSE

c ===  fill /ISOPHO/ with isolated photon 
c ===  remove pho-cluster from /CLUSTER/
			IF(ISOL) {								// usuniecie z przepisaniem
				NPHO = NPHO + 1
				IF(LCLU != 0) {
					DO II from LCLU to NCLU-1 {
						DO III from 1 to 5 {          
							KCLU(II,III) = KCLU(II+1,III) 
							PCLU(II,III) = PCLU(II+1,III)
						}
					}
				NCLU = NCLU - 1
			}

			KPHO(NPHO,1) = NPHO
			KPHO(NPHO,2) = K(I,2)
           		KPHO(NPHO,3) = I
           		KPHO(NPHO,4) = K(K(I,3),2)
           		KPHO(NPHO,5) = 1
           		
			PPHO(NPHO,1) = ETA
           		PPHO(NPHO,2) = PHI
           		PPHO(NPHO,3) = ETA
           		PPHO(NPHO,4) = PHI
           		PPHO(NPHO,5) = PT
		}
       		} // sprawdz czy czegos tu nie zapomniales? wyglada dziwnie
	}

	CALL HF1(IDENT+11,REAL(NPHO),1.0)

c.....arrange photons in falling E_T sequence
	DO IDU from 1 to NPHO {
		ETMAX = 0
		DO IMU from 1 to NPHO {
			IF(KPHO(IMU,5) == 0) continue
			IF(PPHO(IMU,5) < ETMAX) continue
			IMAX = IMU
			ETMAX = PPHO(IMU,5)
		}

		KPHO(IMAX,5) = 0 //used
		DO II from 1 to 5 {
			KDUM(IDU,II) = KPHO(IMAX,II)
			PDUM(IDU,II) = PPHO(IMAX,II)
		}
	}

	DO I from 1 to NPHO {
		DO II from 1 to 5 {
			KPHO(I,II) = KDUM(I,II)
			PPHO(I,II) = PDUM(I,II)
		}
		KPHO(I,1) = I
		KPHO(I,5) = 1
	}

c === check with partons
	IPHO = 0
	IPHISO = 0
	DO I from 1 to NSTOP {
		IF(ABS(K(I,2)) == 22) {
			PT = SQRT(P(I,1)**2+P(I,2)**2) 
			ETA = SIGN(LOG((SQRT(PT**2+P(I,3)**2)+ABS(P(I,3)))/PT),P(I,3)) 
			PHI = ANGLE(P(I,1),P(I,2))
			ENER = 0.0
			ISOL = TRUE
			
			DO II from 1 to NSTOP {
				IF(ABS(K(II,2)) <= 21 && I != II && ABS(K(II,2)) != 12 && ABS(K(II,2)) != 14 && ABS(K(II,2)) != 16) {
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
				IPHO = IPHO + 1

			IF(ABS(ETA) < ETAMAX && PT > PTLMIN && ISOL)
				IPHISO = IPHISO + 1
		}
	}

	CALL HF1(IDENT+21, REAL(IPHO), 1.0)
	CALL HF1(IDENT+31, REAL(IPHISO), 1.0)



C
      ELSEIF(MODE.EQ.1) THEN
C     =========================

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
