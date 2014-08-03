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
	IEVENT++;
/*	
	const vector<Particle>& parts = irecord.particles();
	
	// look for isolated photons, sort cluster commons
	Int32_t NSTOP = 0, NSTART = 1;
	for (int i=0; i<parts.size(); ++i) {
		if (parts[i].stateID != 21) {
			NSTOP = I-1;
			NSTART = I;
			break;
		}
	}

	for (int i=NSTART; i<parts.size(); ++i) {
		if (!parts[i].isStable()) 
			continue;
			
		if (parts[i].pT() == 0) 
			continue;

		// analyse photons
		if (abs(parts[i].typeID) == 22) {
			ISOL = true;
			LCLU = -1;
			PT = parts[i].pT();
			ETA = parts[i].getEta();
			PHI = parts[i].getPhi();
			PTCRU = PT;
			
			// apply smearing
			if (KEYSME) {
				ENE = parts[i].e();
				if (ENE <= 0.0) 
					continue;
					
				SIGPH = RESPHO(ENE,PT,ETA,PHI);
				pxpho = parts[i].pX() * (1.0 + sigph);
				pypho = parts[i].pY() * (1.0 + sigph);
				pzpho = parts[i].pZ() * (1.0 + sigph);
				eepho = parts[i].e()  * (1.0 + sigph);
				PT = sqrt(pxpho*pxpho + pypho*pypho);
				ENE = eepho;
			}

			if (PT < PTLMIN) 
				continue;
				 
			if (abs(ETA) > ETAMAX) 
				continue; 

			// mark photon-cluster
			DR = 100.0;
			for (int j=0; j<orecord.vCluster.size(); ++j) {
				DDR = sqrt(
					pow(ETA - orecord.vCluster[j].P[2], 2) +
					pow(PHI - orecord.vCluster[j].P[3], 2)
				);
				
				if (abs(PHI - orecord.vCluster[j]) > PI)
					DDR = sqrt(
						pow(ETA - orecord.vCluster[j].P[2], 2) + 
						pow(abs(PHI - orecord.vCluster[j].P[3])-2*PI, 2)
					);
					
				if (DDR < DR) 
					LCLU = j;
					
				DR = min(DDR,DR);
			}

			if (DR > RJE) 
				LCLU = -1;
				
			// check if photon  isolated from clusters
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
							pow(abs(PHI - orecord.vCluster[j].P[3])-2*PI, 2)
						);
				}

				if (DDR < RISOLJ) 
					ISOL = false;
			}
		
			// check on energy deposition of cells EDEP in cone RDEP
			Real64_t EDEP = 0.0;
			for (int j=0; j<orecord.vCell.size(); ++j) {
				DDR = sqrt(
					pow(ETA - orecord.vCell[j].P[0], 2) + 
					pow(PHI - orecord.vCell[j].P[1], 2)
				);
				
				if (abs(PHI - orecord.vCell[j].P[1]) > PI)
					DDR = sqrt(
						pow(ETA - orecord.vCell[j].P[0], 2) + 
						pow(abs(PHI - orecord.vCell[j].P[1])-2*PI, 2)
					);
					
				if (DDR < RDEP) 
					EDEP += orecord.vCell[j].P[4];
			}

			if (EDEP - PT > EDMAX) 
				ISOL = false;
		
			// fill /ISOPHO/ with isolated photon 
			// remove pho-cluster from /CLUSTER/
			if (ISOL) {
				if (LCLU != -1) {
					DO II from LCLU to NCLU-1 {
						DO III from 1 to 5 {          
							KCLU(II,III) = KCLU(II+1,III) 
							PCLU(II,III) = PCLU(II+1,III)
						}
					}
					NCLU = NCLU - 1
				}

				KPHO(NPHO,1) = NPHO;
				KPHO(NPHO,2) = K(I,2);
				KPHO(NPHO,3) = I;
				KPHO(NPHO,4) = K(K(I,3),2);
				KPHO(NPHO,5) = 1;
           		
				PPHO(NPHO,1) = ETA;
				PPHO(NPHO,2) = PHI;
				PPHO(NPHO,3) = ETA;
				PPHO(NPHO,4) = PHI;
				PPHO(NPHO,5) = PT;
			}
		}
	}

	// CALL HF1(IDENT+11,REAL(NPHO),1.0)

	// arrange photons in falling E_T sequence
	IDU from 1 to NPHO {
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
	
	// check with partons
	Int32_t IPHO = 0;
	Int32_t IPHISO = 0;
	for (int i=0; i<=NSTOP; ++i) {
		if (abs(parts[i].typeID) == 22) {
			PT = parts[i].pT(); 
			ETA = parts[i].getEta();
			PHI = parts[i].getPhi();
			ENER = 0.0;
			ISOL = true;
			
			for (int j=0; j<=NSTOP; ++j) {
				if (abs(parts[j].typeID) <= 21 && i != j && 
					abs(parts[j].typeID) != 12 && 
					abs(parts[j].typeID) != 14 && 
					abs(parts[j].typeID) != 16) 
				{
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
							pow(abs(JPHI-PHI)-2*PI, 2)
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
				IPHO++;

			if (abs(ETA) < ETAMAX && PT > PTLMIN && ISOL)
				IPHISO++;
		}
	}

	// CALL HF1(IDENT+21, REAL(IPHO), 1.0)
	// CALL HF1(IDENT+31, REAL(IPHISO), 1.0)
*/
}

void Photon::printResults() const {
	printf ("*************************************\n");
	printf ("*                                   *\n");
	printf ("*     *************************     *\n");
	printf ("*     ***    Output from    ***     *\n");
	printf ("*     ***  analyse::Photon  ***     *\n");
	printf ("*     *************************     *\n");
	printf ("*                                   *\n");
	printf ("*************************************\n");
	
	printf (" Analysed records: %d\n", IEVENT);
	histo_isol	.print( true );
	histo_hard	.print( true );
	histo_sum	.print( true );
}
