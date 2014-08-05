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
				const ClusterData& cluster = orecord.Clusters[j];
				
				DDR = sqrt(
					pow(ETA - cluster.eta_rec, 2) +
					pow(PHI - cluster.phi_rec, 2)
				);
				
				if (abs(PHI - cluster.phi_rec) > PI)
					DDR = sqrt( 
						pow(ETA - cluster.eta_rec, 2) + 
						pow(abs(PHI - cluster.phi_rec) - 2*PI), 2)
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
					const ClusterData& cluster = orecord.Clusters[j];
					
					DDR = sqrt(
						pow(ETA - cluster.eta_rec, 2) + 
						pow(PHI - cluster.phi_rec, 2)
					);
					
					if (abs(PHI - cluster.phi_rec) > PI)
						DDR = sqrt(
							pow(ETA - cluster.eta_rec, 2) + 
							pow(abs(PHI - cluster.phi_rec) - 2*PI, 2)
						);
				}

				if (DDR < RISOLJ) 
					ISOL = false;
			}
			
			// check on energy deposition of cells EDEP in cone RDEP
			EDEP = 0.0;
			for (int j=0; j<orecord.Cells.size(); ++j) {
				const CellData& cell = orecord.Cells[j];
				
				DDR = sqrt(
					pow(ETA - cell.eta, 2) + 
					pow(PHI - cell.phi, 2)
				);
				
				if (abs(PHI - cell.phi) > PI)
					DDR = sqrt(
						pow(ETA - cell.eta, 2) + 
						pow(abs(PHI - cell.phi) - 2*PI, 2)
					);
					
				if (DDR < RDEP) 
					EDEP += cell.pT;
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
				KELE(NELE,2) = K(I,2)			// typ elektronu
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
	// sor KELE i PELE

	Int32_t IELE = 0, IELEISO = 0;
	Real64_t PT, ETA, PHI, ENER, JPT, JETA, JPHI, DDR;
	for (int i=0; i<=NSTOP; ++i) {
		const Particle& part = parts[i];
		
		if (part.type == PT_ELECTRON) {
			PT = part.pT(); 
			ETA = part.getEta(); 
			PHI = part.getPhi();
			ENER = 0.0;
			Bool_t ISOL = true;

			for (int j=0; j<=NSTOP; ++j) {
				if (abs(parts[j].typeID) <= 21 && i != j && !parts[j].isNeutrino()) {
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
