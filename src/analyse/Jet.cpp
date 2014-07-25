#include "Jet.h"
#include <cstdio>

using namespace AcerDet::analyse;

Jet::Jet( const Configuration& config ) :
	ETJET	( config.Jet.MinEnergy ),
	ETAJET	( config.Jet.RapidityCoverage ),
	RCONE	( config.Cluster.ConeR ),
	PTMIN	( config.Cell.MinpT ),
	CALOTH	( config.Cell.EtaTransition ),

	KEYHID	( config.Flag.HistogramId ),
	KEYSME	( config.Flag.Smearing ),
	KEYFLD	( config.Flag.BField ),
	KFINVS	( config.Flag.SusyParticle ),

	IEVENT	( 0 ),
	
	histo_bJets				("Jet: multiplicity", 0.0, 10.0, 10),
	histo_delta_phi			("Jet: delta phi jet-barycentre", -0.5, 0.5, 50),
	histo_delta_eta			("Jet: delta eta jet-barycentre", -0.5, 0.5, 50),
	histo_delta_barycenter	("Jet: delta r jet-barycentre", 0.0, 0.5, 50),
	histo_delta_parton		("Jet: delta r jet-parton", 0.0, 0.5, 50),
	histo_pT_bySum			("Jet: pTjet / SumpTParticle", 0.0, 2.0, 50),
	histo_pT_byPart			("Jet: pTjet / pTparton", 0.0, 2.0, 50)
{}

Jet::~Jet() {}

void Jet::printInfo() const {
	// print out title
	printf ("************************************\n");
	printf ("*                                  *\n");
	printf ("*     ************************     *\n");
	printf ("*     ***   analyse::Jet   ***     *\n");
	printf ("*     ************************     *\n");
	printf ("*                                  *\n");
	printf ("************************************\n");

	// print out basic params
	printf ("\n\t clusters definition ....\n");
	printf (" R cone %lf\n", RCONE);
	printf ("\t jets definition ....\n");
	printf (" E_T_jets [GeV] %lf\n", ETJET);
	printf (" eta coverage jets %lf\n", ETAJET);
    printf (" smearing %s\n", KEYSME ? "on" : "off");
    printf (" B-field %s\n", KEYFLD ? "on" : "off");
	printf ("\n");
}

void Jet::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	IEVENT++;
	
	// smear clusters energy
	if (KEYSME) {      
		for (int i=0; i<orecord.vCluster.size(); ++i) {
			EECLU = orecord.vCluster[i].P[4] * cosh(orecord.vCluster[i].P[2]);
			SIGMA = RESHAD(EECLU,orecord.vCluster[i].P[2],CALOTH,orecord.vCluster[i].P[4],RCONE);
			
			Real64_t coef = orecord.vCluster[i].P[4] * (1.0 + SIGMA);
			PXCLU = coef * cos (orecord.vCluster[i].P[3]);
			PYCLU = coef * sin (orecord.vCluster[i].P[3]);
			PZCLU = coef * sinh(orecord.vCluster[i].P[2]);
			EECLU = coef * cosh(orecord.vCluster[i].P[2]);

			// PCLU.getTheta()   az sie prosi zapisac wektorowo
			if (abs(PZCLU / EECLU) < 1) {
				THETA = ACOS(PZCLU / SQRT(PXCLU**2+PYCLU**2+PZCLU**2))
			} ELSE {
				IF (PZCLU > 0) THETA = 0
				IF (PZCLU < 0) THETA = PI
			}

			ETACLU = -log(max(0.0001, abs(tan(0.5*THETA))));
			PTCLU  = SQRT(PXCLU**2+PYCLU**2)  // pClu.pT()
			PHICLU = ANGLE(PXCLU,PYCLU) // pClu.getPhi();

			PCLU(I,5) = PTCLU
			PCLU(I,3) = ETACLU
			PCLU(I,4) = PHICLU
		}
	}
	
	/*
c.....add nonisolated muons to jets
	DO I from 1 to NMUOX {
		ETA = PMUOX(I,3)
		PHI = PMUOX(I,4)
		PT = PMUOX(I,5)

		DR = 100.0
		MUCLU = 0
		DO II from 1 to NCLU {
			DDR = SQRT((ETA-PCLU(II,3))**2+(PHI-PCLU(II,4))**2)
			IF(ABS(PHI-PCLU(II,4)) > PI)
				DDR = SQRT( (ETA-PCLU(II,3))**2 + (ABS(PHI-PCLU(II,4))-2*PI)**2 )

			IF(DDR < DR) {
				MUCLU = II
				DR = DDR
			}
		}

		IF(MUCLU != 0) {
			IF((ABS(PCLU(MUCLU,3)) <= CALOTH && DR < RCONE) || (ABS(PCLU(MUCLU,3)) > CALOTH && DR < RCONE) ) {
				PXCLU = PCLU(MUCLU,5)*COS(PCLU(MUCLU,4)) + PT*COS(PHI)
				PYCLU = PCLU(MUCLU,5)*SIN(PCLU(MUCLU,4)) + PT*SIN(PHI)
				PZCLU = PCLU(MUCLU,5)*SINH(PCLU(MUCLU,3)) + PT*SINH(ETA)
				EECLU = PCLU(MUCLU,5)*COSH(PCLU(MUCLU,3)) + PT*COSH(ETA)

				PTCLU = SQRT(PXCLU**2+PYCLU**2)
				ETACLU = SIGN(LOG((SQRT(PTCLU**2+PZCLU**2)+ABS(PZCLU))/PTCLU),PZCLU) 
				PHICLU = ANGLE(PXCLU,PYCLU)

				PCLU(MUCLU,3) = ETACLU
				PCLU(MUCLU,4) = PHICLU
				PCLU(MUCLU,5) = PTCLU

				KMUOX(I,5) = 0
			}
		}
	}

c.... store accepted jets in \JETALL\ common and flagg in common /CLUSTER/
	NJET = 0
	DO I from 1 to NCLU {
		IF(PCLU(I,5) > ETJET && ABS(PCLU(I,3)) < ETAJET) {
			NJET = NJET + 1
			DO II from 1 to 5 {          
				KJET(NJET,II) = KCLU(I,II) 
				PJET(NJET,II) = PCLU(I,II)
			}

			KJET(NJET,1) = NJET
			KCLU(I,5) = 0
		}
	}

c.....histogram NJET
	CALL HF1(IDENT+1,REAL(NJET),1.0)

c.....arrange jets in falling E_T sequence // sort by E_T
	DO IDU from 1 to NJET {
		ETMAX = 0
		DO IMU from 1 to NJET {
			IF(KJET(IMU,5) == 0) continue
			IF(PJET(IMU,5) < ETMAX) continue
			IMAX = IMU
			ETMAX = PJET(IMU,5)
		}

		KJET(IMAX,5) = 0
		DO II from 1 to 5 {
			KDUM(IDU,II) = KJET(IMAX,II)
			PDUM(IDU,II) = PJET(IMAX,II)
		}
	}

	DO I from 1 to NJET {
		DO II from 1 to 5 {
			KJET(I,II) = KDUM(I,II)
			PJET(I,II) = PDUM(I,II)
		}
		KJET(I,1) = I
		KJET(I,5) = 1
	}

c....reconstruct baricenter of particles
	DO IJET from 1 to NJET {
		ETAREC = 0
		PTREC = 0
		PHIREC = 0

		DO I from 1 to N { 
			IF(K(I,1) <= 0 || K(I,1) > 10) continue
			IF(P(I,1)**2+P(I,2)**2 <= PTLRAT*P(I,3)**2) continue 

			KC = KFCOMP(K(I,2)) 
			IF(KC == 0 || KC == 12 || KC == 14 || KC == 16 || KC == 13 || KC == KFINVS) continue 

			DETPHI = 0.0
			IF(KEYFLD && KUCHGE(K(I,2)) != 0) {
				PT = SQRT(P(I,1)**2+P(I,2)**2)
				IF(PT < PTMIN) 
					continue

				ETA = SIGN(LOG((SQRT(PT**2+P(I,3)**2)+ABS(P(I,3)))/PT),P(I,3)) 
				PHI = ANGLE(P(I,1),P(I,2))
				CHRG = KUCHGE(K(I,2))/3.

				IF(ABS(P(I,3)/P(I,4)) < 1) {
					THETA = ACOS(P(I,3)/SQRT(PT**2+P(I,3)**2))
				} ELSE {
					IF(P(I,3) > 0) THETA = 0
					IF(P(I,3) < 0) THETA = PI
				}

				DETPHI = CHRG*FLDPHI(PT,ETA,THETA)
			}

			PT = SQRT(P(I,1)**2+P(I,2)**2)
			ETA = SIGN(LOG((SQRT(PT**2+P(I,3)**2)+ABS(P(I,3)))/PT),P(I,3)) 
			PHI = ANGLE(P(I,1),P(I,2))
			PHI = PHI + DETPHI
			DPHIA = ABS(PJET(IJET,2)-PHI) 

			IF(DPHIA > PI) 
				DPHIA = DPHIA - 2*PI

			IF(ABS(PJET(IJET,1)) < CALOTH && (PJET(IJET,1)-ETA)**2+(DPHIA)**2 > RCONE**2) continue 
			IF(ABS(PJET(IJET,1)) > CALOTH && (PJET(IJET,1)-ETA)**2+(DPHIA)**2 > RCONE**2) continue 

			PTREC = PTREC + PT
			ETAREC = ETAREC + ETA * PT
			PHIREC = PHIREC + PHI * PT
		}

		ETAREC /= PTREC
		PHIREC /= PTREC
		DETR = SQRT((ETAREC-PJET(IJET,3))**2+(PHIREC-PJET(IJET,4))**2)

		CALL HF1(IDENT+11,ETAREC-PJET(IJET,3),1.0)
		CALL HF1(IDENT+12,PHIREC-PJET(IJET,4),1.0)
		CALL HF1(IDENT+13,DETR,1.0)
		CALL HF1(IDENT+14,PJET(IJET,5)/PTREC,1.0)
	}

	DO IJET from 1 to NJET {
		PTREC = 0
		DETRMIN = RCONE

		DO I from 7 to N { 
			IF(K(I,1) != 21 || ABS(K(I,2)) > 10) continue

			PT = SQRT(P(I,1)**2+P(I,2)**2)
			ETA = SIGN(LOG((SQRT(PT**2+P(I,3)**2)+ABS(P(I,3)))/PT),P(I,3)) 
			PHI = ANGLE(P(I,1),P(I,2))

			DPHIA = ABS(PJET(IJET,4)-PHI) 
			IF(DPHIA > PI) 
				DPHIA = DPHIA - 2*PI

			DETR = SQRT((ETA-PJET(IJET,3))**2+DPHIA**2)
			IF(DETR > DETRMIN) continue
			
			PTREC = PT
			DETRMIN = DETR
		}

		IF(PTREC != 0) {
			CALL HF1(IDENT+23,DETRMIN,1.0)
			CALL HF1(IDENT+24,PJET(IJET,5)/PTREC,1.0)
		}
	}


C
      ELSEIF(MODE.EQ.1) THEN
C     =========================

      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '        OUTPUT FROM AcerDET      '
      WRITE(NOUT,BXTXT) '              ACDJET             '
      WRITE(NOUT,BXTXT) '*********************************'
      
      WRITE(NOUT,BXCLO)

      CALL HPRINT(IDENT+1)
      CALL HPRINT(IDENT+11)
      CALL HPRINT(IDENT+12)
      CALL HPRINT(IDENT+13)
      CALL HPRINT(IDENT+14)
      CALL HPRINT(IDENT+23)
      CALL HPRINT(IDENT+24)

      ENDIF
C     =====

      END
*/

}
