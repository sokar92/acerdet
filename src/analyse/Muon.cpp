#include "Muon.h"
#include <cstdio>

using namespace AcerDet::analyse;

Muon::Muon( const Configuration& config ) :
	ETCLU	( config.Cluster.MinEt ),
	RCONE	( config.Cluster.ConeR ),

	PTMUMIN	( config.Muon.MinMomenta ),
	ETAMAX	( config.Muon.MaxEta ),
	RISOLJ	( config.Muon.MinIsolRlj ),
	RDEP	( config.Muon.ConeR ),
	EDMAX	( config.Muon.MaxEnergy ),
	PTMUMINT( config.Muon.MinMomenta ),    // PROBLEM!!!

	KEYHID	( config.Flag.HistogramId ),
	KEYSME	( config.Flag.Smearing ),

	IEVENT	( 0 ),
	
	histo_nonisol	("Muon: non-isolated", 0.0, 10.0, 10),
	histo_isol		("Muon: isolated", 0.0, 10.0, 10),
	histo_hard		("Muon: hard", 0.0, 10.0, 10),
	histo_sum		("Muon: hard+isol", 0.0, 10.0, 10)
{}

Muon::~Muon() {}

void Muon::printInfo() const {
	// print out title
	printf ("*************************************\n");
	printf ("*                                   *\n");
	printf ("*     *************************     *\n");
	printf ("*     ***   analyse::Muon   ***     *\n");
	printf ("*     *************************     *\n");
	printf ("*                                   *\n");
	printf ("*************************************\n");

	// print out basic params
	printf ("\n\t... muon isolation ...\n");
	printf ("min. muon p_T %lf\n", PTMUMIN);
	printf ("max. muon eta %lf\n", ETAMAX);
	printf ("min R_lj for isolation %lf\n", RISOLJ);
	printf ("R for energy deposit %lf\n", RDEP);
	printf ("max E_dep for isolation %lf\n", EDMAX);
	printf ("min muon p_T unsmea %lf\n", PTMUMINT);
	printf ("smearing %s\n", KEYSME ? "on" : "off");
	printf ("\n");
}

void Muon::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	//printf ("Muon: analyse record\n");
	/*
      SUBROUTINE ACDMUO(MODE,LPAR,YPAR)

      INTEGER  KEYHID, KEYSME, IEVENT 
      INTEGER I, II, NSTART, NSTOP, IMUO, IMUOISO
      REAL DDR, EDEP
      REAL ANGLE
      REAL ETA, PHI, PT
      REAL RISOLJ, RDEP, EDMAX, PTMUMIN, ETCLU, PTMUMINT
      REAL RCONE, ETAMAX
      REAL RESMUO, SIGMA
      REAL PXMUO, PYMUO, PZMUO, EEMUO
      REAL ENER, JPT, JPHI, JETA
      LOGICAL ISOL
      INTEGER NBINA, IDENT
      REAl TMAXA,TMINA
      REAL YPAR
      INTEGER MODE, LPAR 

      ELSEIF(MODE.EQ.0) THEN
C     ======================
	IEVENT = IEVENT+1

	// wyczysc tablice przed przetwarzaniem (nic nie rob)
	NMUO = 0
	DO I=1,100 {
		KMUO(I,2)=0
	}

	NMUOX = 0
	DO I=1,100 {
		KMUOX(I,2) = 0
	}

	// znajdz poczatek danych
	NSTOP = 0
	NSTART = 1
	DO I=1, N {
		IF(K(I,1) != 21) 
		THEN
			NSTOP = I-1
			NSTART= I
			break
		ENDIF
	}

c....look for isolated muons, sort clusters common
	DO I from NSTART to N {
		IF(K(I,1) < 0 || K(I,1) > 10) continue
		IF(ABS(K(I,2)) == 13) 					// czastka jest muonem
		THEN

c -- analyse not decayed muons
			ISOL = TRUE
			PT = SQRT(P(I,1)**2+P(I,2)**2)
			IF(PT < PTMUMINT) continue
			ETA = SIGN(LOG((SQRT(PT**2+P(I,3)**2)+ABS(P(I,3)))/PT),P(I,3)) 
			IF(ABS(ETA) > ETAMAX) continue  
			PHI = ANGLE(P(I,1),P(I,2))

c.... apply smearing
			IF(KEYSME) 
			THEN
				SIGMA = RESMUO(PT,ETA,PHI)
				PXMUO = P(I,1) / (1.+SIGMA)
				PYMUO = P(I,2) / (1.+SIGMA)
				PZMUO = P(I,3) / (1.+SIGMA)
				EEMUO = P(I,4) / (1.+SIGMA)
				PT    = SQRT(PXMUO**2+PYMUO**2)
			ENDIF
			IF(PT < PTMUMIN) continue

c ===   check is muon isolated from clusters
			DO II from 1 to NCLU {
				DDR = SQRT((ETA-PCLU(II,3))**2+(PHI-PCLU(II,4))**2)
				IF(ABS(PHI-PCLU(II,4)) > PI)
					DDR = SQRT( (ETA-PCLU(II,3))**2 + (ABS(PHI-PCLU(II,4))-2*PI)**2 )
				IF(DDR < RISOLJ) 
					ISOL = FALSE
			}

c ===   check on energy deposition of   cells EDEP in cone RDEP
			EDEP = 0.0
			DO II from 1 to NCELL {
				DDR = SQRT((ETA-PCELL(II,1))**2+(PHI-PCELL(II,2))**2)
				IF(ABS(PHI-PCELL(II,2)) > PI)
					DDR=SQRT( (ETA-PCELL(II,1))**2 + (ABS(PHI-PCELL(II,2))-2*PI)**2 )
				IF(DDR < RDEP) 
					EDEP = EDEP + PCELL(II,5)
			}
			IF(EDEP > EDMAX) 
				ISOL = FALSE

c ===   fill /ISOMUO/ with isolated muon
			IF(ISOL) 
			THEN
				NMUO = NMUO + 1
				KMUO(NMUO,1) = NMUO
				KMUO(NMUO,2) = K(I,2)
				KMUO(NMUO,3) = I
				KMUO(NMUO,4) = K(K(I,3),2)
				KMUO(NMUO,5) = 1

				PMUO(NMUO,1) = ETA
				PMUO(NMUO,2) = PHI
				PMUO(NMUO,3) = ETA
				PMUO(NMUO,4) = PHI
				PMUO(NMUO,5) = PT
			ELSE
c ===   fill /NOISOMUO/ with non-isolated muon
				NMUOX = NMUOX + 1
				KMUOX(NMUOX,1) = NMUOX
				KMUOX(NMUOX,2) = K(I,2)
          			KMUOX(NMUOX,3) = I
				KMUOX(NMUOX,4) = K(K(I,3),2)
				KMUOX(NMUOX,5) = 1

          			PMUOX(NMUOX,1) = ETA
          			PMUOX(NMUOX,2) = PHI
          			PMUOX(NMUOX,3) = ETA
          			PMUOX(NMUOX,4) = PHI
          			PMUOX(NMUOX,5) = PT
			ENDIF //isol
		ENDIF //is muon
c.............
	} // end for

	CALL HF1(IDENT+11,REAl(NMUO),1.0)
	CALL HF1(IDENT+10,REAl(NMUOX),1.0)

c....cross-check with partons
	IMUO = 0
	IMUOISO = 0
	DO I from 1 to NSTOP {
		IF(ABS(K(I,2)) == 13) {
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
					DDR = SQRT((ETA-JETA)**2 + (JPHI-PHI)**2)
					IF(ABS(JPHI-PHI) > PI)
						DDR = SQRT( (ETA-JETA)**2 + (ABS(JPHI-PHI)-2*PI)**2 )
					IF(DDR < RISOLJ && JPT > ETCLU) 
						ISOL = FALSE
					IF(DDR < RDEP && JPT < ETCLU) 
						ENER = ENER + JPT
				}
			}

			IF(ENER > EDMAX) 
				ISOL = FALSE
			IF(ABS(ETA) < ETAMAX && PT > PTMUMIN) 
				IMUO = IMUO + 1
			IF(ABS(ETA) < ETAMAX && PT > PTMUMIN && ISOL)
				IMUOISO = IMUOISO + 1 
		}
	}

	CALL HF1(IDENT+21,REAL(IMUO),1.0)
	CALL HF1(IDENT+31,REAL(IMUOISO),1.0)

C
      ELSEIF(MODE.EQ.1) THEN
C     =========================

      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '        OUTPUT FROM AcerDET      '
      WRITE(NOUT,BXTXT) '              ACDMUO             '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXCLO)

      CALL HPRINT(IDENT+10)
      CALL HPRINT(IDENT+11)
      CALL HPRINT(IDENT+21)
      CALL HPRINT(IDENT+31)

      ENDIF
C     =====

      END
*/

}
