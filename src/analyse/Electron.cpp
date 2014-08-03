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
	IEVENT++;
/*	
	const vector<Particle>& parts;
	
	Int32_t NSTOP = 0, NSTART = 1;
	for (int i=0; i<parts.size(); ++i) {
		if (parts[i].statusID != 21) {
			NSTOP = i-1;
			NSTART = i;
			break;
		}
	}
	
	// look for isolated electrons, sort clusters common
	for (int i=NSTART; i<parts.size(); ++i) {
		if (!parts[i].isStable()) 
			continue;
			
		// analyse electrons
		if (abs(parts[i].typeID) == 11) {
			ISOL = true;
			LCLU = 0;
			PT = parts[i].pT();
        	ETA = parts[i].getEta();
        	PHI = parts[i].getPhi();
   
			// apply smearing
        	PTCRU = PT;
			if (KEYSME) {
				ENE = parts[i].e();
				if (ENE <= 0.0) 
					continue;

				SIGPH = RESELE(ENE,PT,ETA,PHI);
				PXELE = P(I,1) * (1.0 + SIGPH);
				PYELE = P(I,2) * (1.0 + SIGPH);
				PZELE = P(I,3) * (1.0 + SIGPH);
				EEELE = P(I,4) * (1.0 + SIGPH);
				PT = sqrt(PXELE*PXELE + PYELE*PYELE);
			}
        
			if (PT < PTLMIN)
				continue;

			if (abs(ETA) > ETAMAX)
				continue;
				
			// mark eletron-cluster
			DR = 100.0;
			for (int j=0; j<orecord.vCluster.size(); ++j) {
				DDR = sqrt(
					pow(ETA - orecord.vCluster[j].P[2], 2) +
					pow(PHI - orecord.vCluster[j].P[3], 2)
				);
				
				if (abs(PHI - orecord.vCluster[j].P[3]) > PI)
					DDR = sqrt( 
						pow(ETA - orecord.vCluster[j].P[2], 2) + 
						pow(abs(PHI - orecord.vCluster[j].P[3]) - 2*PI), 2)
					);
					
				if (DDR < DR) 
					LCLU = j;
					
				DR = MIN(DDR,DR);
			}

			if (DR > RJE) 
				LCLU = -1;

			// check if electron isolated from clusters
			for (int j=0; j<orecord.vCluster.size(); ++j) {
				if (j == LCLU) {
					DDR = 100.0;
				} else {
					DDR = sqrt(
						pow(ETA - orecord.vCluster[j].P[2], 2) + 
						pow(PHI - orecord.vCluster[j].P[3], 2)
					);
					
					if (abs(PHI - orecord.vCluster[j].P[3]) > PI)
						DDR = sqrt(
							pow(ETA - orecord.vCluster[j].P[2], 2) + 
							pow(abs(PHI - orecord.vCluster.P[3]) - 2*PI, 2)
						);
				}

				if (DDR < RISOLJ) 
					ISOL = false;
			}
			
			// check on energy deposition of cells EDEP in cone RDEP
			EDEP = 0.0;
			for (int j=0; j<orecord.vCell.size(); ++j) {
				DDR = sqrt(
					pow(ETA - orecord.vCell[j].P[0], 2) + 
					pow(PHI - orecord.vCell[j].P[1], 2)
				);
				
				if (abs(PHI - orecord.vCell[j].P[1]) > PI)
					DDR = sqrt(
						pow(ETA - orecord.vCell[j].P[0], 2) + 
						pow(abs(PHI - orecord.vCell[j].P[1]) - 2*PI, 2)
					);
					
				if (DDR < RDEP) 
					EDEP += orecord.vCell[j].P[4];
			}

			if (EDEP - PTCRU > EDMAX) 
				ISOL = false;
				
			// fill /ISOELE/ with isolated electron 
			// remove ele-cluster from /CLUSTER/
			if (ISOL) {
				 usun stowarzyszony cluster z listy clustrow
				if (LCLU != -1) {
					DO II from LCLU to NCLU-1 {
						DO III from 1 to 5 {
							KCLU(II,III) = KCLU(II+1,III) 
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

	// CALL HF1(IDENT+11, REAL(NELE), 1.0)

	// arrange electrons in falling E_T sequence
	DO IDU from 1 to NELE {
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
	Int32_t IELE = 0;
	Int32_t IELEISO = 0;

	for (int i=0; i<=NSTOP; ++i) {
		if (abs (parts[i].typeID) == 11) {
			PT = parts[i].pT(); 
			ETA = parts[i].getEta(); 
			PHI = parts[i].getPhi();
			ENER = 0.0;
			ISOL = true;

			for (int j=0; j<=NSTOP; ++j) {
				if (abs(parts[j].typeID) <= 21 && i != j && 
					abs(parts[j].partID) != 12 && 
					abs(parts[j].partID) != 14 && 
					abs(parts[j].partID) != 16) 
				{
					// przelicz dane dla II czastki - standard per particle
					JPT = parts[j].pT(); 
					JETA = parts[j].getEta(); 
					JPHI = parts[j].getPhi();
					DDR = sqrt(
						pow(ETA-JETA, 2) + 
						pow(JPHI-PHI, 2)
					);

					if (abs(JPHI-PHI) > PI)
						DDR = sqrt(
							pow(ETA-JETA, 2) + 
							pow(abs(JPHI-PHI)-2*PI), 2)
						);

					if (DDR < RISOLJ && JPT > ETCLU) 
						ISOL = false;

					if (DDR < RDEP && JPT < ETCLU) 
						ENER += JPT;
				}
			}

			if (ENER > EDMAX) 
				ISOL = false;

			if (abs(ETA) < ETAMAX && PT > PTLMIN)
				IELE++;

			if (abs(ETA) < ETAMAX && PT > PTLMIN && ISOL)
				IELEISO++;
		}
	}
	
	// CALL HF1(IDENT+21, REAL(IELE), 1.0)
	// CALL HF1(IDENT+31, REAL(IELEISO), 1.0)
*/
}

void Electron::printResults() const {
	printf ("***************************************\n");
	printf ("*                                     *\n");
	printf ("*     ***************************     *\n");
	printf ("*     ***     Output from     ***     *\n");
	printf ("*     ***  analyse::Electron  ***     *\n");
	printf ("*     ***************************     *\n");
	printf ("*                                     *\n");
	printf ("***************************************\n");
	
	printf (" Analysed records: %d\n", IEVENT);
	histo_isol	.print( true );
	histo_hard	.print( true );
	histo_sum	.print( true );
}
