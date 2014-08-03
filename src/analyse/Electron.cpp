#include "Electron.h"
#include <cstdio>

using namespace AcerDet::analyse;

Electron::Electron( const Configuration& config ) :
	ETCLU	( config.Cluster.MinEt ),
	RCONE	( config.Cluster.ConeR ),

	PTLMIN	( config.Electron.MinMomenta ),
	ETAMAX	( config.Electron.MaxEta ),
	RJE		( config.Electron.MinJetsRlj ),
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
	/*
	// new event to compute
	IEVENT++;

	// reference to particles container
	const vector<Particle>& parts = irecord.particles();

	// znajdz poczatek danych
	Int32_t NSTOP = 0, NSTART = 1;
	for (int i=0; i<parts.size(); ++i) {
		if (parts[i].stateID != 21) {
			NSTOP = i-1;
			NSTART = i;
			break;
		}
	}

	// look for isolated electrons, sort clusters common
	for (int i=NSTART; i<parts.size(); ++i) {
		const Particle& part = parts[i];
	
		if (!part.isStable()) 
			continue;
			
		// analyse electrons
		if (part.type == PT_ELECTRON) {
			Bool_t ISOL = true;
			LCLU = 0;
			
			PT = part.pT();
        	ETA = part.getEta();
        	PHI = part.getPhi();
   
			// apply smearing
        	PTCRU = PT;
			if (KEYSME) {
				ENE = part.e();
				if (ENE <= 0.0) 
					continue;

				// SIGPH = RESELE(ENE,PT,ETA,PHI);
				// PXELE = P(I,1) * (1.0 + SIGPH);
				// PYELE = P(I,2) * (1.0 + SIGPH);
				// PZELE = P(I,3) * (1.0 + SIGPH);
				// EEELE = P(I,4) * (1.0 + SIGPH);
				// PT = sqrt(PXELE*PXELE + PYELE*PYELE);
				
				Particle pEle = part;
				pEle.momentum *= (1.0 + SIGPH);
				PT = pEle.pT();
			}
        
			if (PT < PTLMIN)
				continue;

			if (abs(ETA) > ETAMAX)
				continue;
				
			// mark eletron-cluster
			DR = 100.0;
			for (int j=0; j<orecord.Clusters.size(); ++j) {
				DDR = sqrt(
					pow(ETA - orecord.Clusters[j].eta_rec, 2) +
					pow(PHI - orecord.Clusters[j].phi_rec, 2)
				);
				
				if (abs(PHI - orecord.Clusters[j].phi_rec) > PI)
					DDR = sqrt( 
						pow(ETA - orecord.Clusters[j].eta_rec, 2) + 
						pow(abs(PHI - orecord.Clusters[j].phi_rec) - 2*PI), 2)
					);
					
				if (DDR < DR) {
					LCLU = j;
					DR = DDR;
				}
			}

			if (DR > RJE) 
				LCLU = -1;

			// check if electron isolated from clusters
			for (int j=0; j<orecord.Clusters.size(); ++j) {
				if (j == LCLU) {
					DDR = 100.0;
				} else {
					DDR = sqrt(
						pow(ETA - orecord.Clusters[j].eta_rec, 2) + 
						pow(PHI - orecord.Clusters[j].phi_rec, 2)
					);
					
					if (abs(PHI - orecord.Clusters[j].phi_rec) > PI)
						DDR = sqrt(
							pow(ETA - orecord.Clusters[j].eta_rec, 2) + 
							pow(abs(PHI - orecord.Clusters[j].phi_rec) - 2*PI, 2)
						);
				}

				if (DDR < RISOLJ) 
					ISOL = false;
			}
			
			// check on energy deposition of cells EDEP in cone RDEP
			EDEP = 0.0;
			for (int j=0; j<orecord.Cells.size(); ++j) {
				DDR = sqrt(
					pow(ETA - orecord.Cells[j].eta, 2) + 
					pow(PHI - orecord.Cells[j].phi, 2)
				);
				
				if (abs(PHI - orecord.Cells[j].phi) > PI)
					DDR = sqrt(
						pow(ETA - orecord.Cells[j].eta, 2) + 
						pow(abs(PHI - orecord.Cells[j].phi) - 2*PI, 2)
					);
					
				if (DDR < RDEP) 
					EDEP += orecord.Cells[j].pT;
			}

			if (EDEP - PTCRU > EDMAX) 
				ISOL = false;
				
			// fill /ISOELE/ with isolated electron 
			// remove ele-cluster from /CLUSTER/
			if (ISOL) {
				// usun stowarzyszony cluster z listy clustrow
				if (LCLU >= 0)
					orecord.Cluster.erase(orecord.Cluster.begin() + LCLU);

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
	//DO IDU from 1 to NELE {
	//	ETMAX = 0
	//	DO IMU from 1 to NELE {
	//		IF(KELE(IMU,5) == 0) continu
	//		IF(PELE(IMU,5) < ETMAX) continue
	//		IMAX = IMU
	//		ETMAX = PELE(IMU,5)
	//	}

	//	KELE(IMAX,5) = 0
	//	DO II from 1 to 5 {
	//		KDUM(IDU,II) = KELE(IMAX,II)
	//		PDUM(IDU,II) = PELE(IMAX,II)
	//	}
	//}

	//DO I from 1 to NELE {
	//	DO II from 1 to 5 {				// skopiuj dum -> ele
	//		KELE(I,II) = KDUM(I,II)
	//		PELE(I,II) = PDUM(I,II)
	//	}
	//	KELE(I,1) = I					// numer elektrona
	//	KELE(I,5) = 1					// stan(1) = przetworzony?
	//}

	Int32_t IELE = 0, IELEISO = 0;
	for (int i=0; i<=NSTOP; ++i) {
		const Particle& part = parts[i];
		
		if (part.type == PT_ELECTRON) {
			PT = part.pT(); 
			ETA = part.getEta(); 
			PHI = part.getPhi();
			ENER = 0.0;
			ISOL = true;

			for (int j=0; j<=NSTOP; ++j) {
				if (abs(parts[j].typeID) <= 21 && i != j && 
					abs(parts[j].partID) != 12 &&  // TODO zamien to na typy enumowe -> zwieksz czytelnosc 
					abs(parts[j].partID) != 14 && 
					abs(parts[j].partID) != 16) 
				{
					JPT = parts[j].pT(); 
					JETA = parts[j].getEta(); 
					JPHI = parts[j].getPhi();
					DDR = sqrt(
						pow(ETA - JETA, 2) + 
						pow(JPHI - PHI, 2)
					);

					if (abs(JPHI - PHI) > PI)
						DDR = sqrt(
							pow(ETA - JETA, 2) + 
							pow(abs(JPHI - PHI) - 2*PI), 2)
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
