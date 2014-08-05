#include "Photon.h"
#include <cstdio>

using namespace AcerDet::analyse;

Photon::Photon( const Configuration& config ) :
	ETCLU	( config.Cluster.MinEt ),
	RCONE	( config.Cluster.ConeR ),

	PTLMIN	( config.Photon.MinMomenta ),
	ETAMAX	( config.Photon.MaxEta ),
	RJE		( config.Photon.MinJetsRlj ),
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
		
		if (parts.pT() == 0) 
			continue;
			
		// analyse photons
		if (part.type == PT_PHOTON) {
			Boolt_t ISOL = true;
			Int32_t LCLU = -1;
			
			PT = part.pT();
			ETA = part.getEta();
			PHI = part.getPhi();
			PTCRU = PT;
			
			// apply smearing
			if (KEYSME) {
				ENE = parts[i].e();
				if (ENE <= 0.0) 
					continue;
					
				// SIGPH = RESPHO(ENE,PT,ETA,PHI);
				// pxpho = parts[i].pX() * (1.0 + sigph);
				// pypho = parts[i].pY() * (1.0 + sigph);
				// pzpho = parts[i].pZ() * (1.0 + sigph);
				// eepho = parts[i].e()  * (1.0 + sigph);
				// PT = sqrt(pxpho*pxpho + pypho*pypho);
				// ENE = eepho;
				
				Particle pPho = part;
				pPho.momentum *= (1.0 + SIGPH);
				
				PT = pPho.pT();
				ENE = pPho.e();
			}

			if (PT < PTLMIN) 
				continue;
				 
			if (abs(ETA) > ETAMAX) 
				continue; 

			// mark photon-cluster
			Real64_t DR = 100.0;
			for (int j=0; j<orecord.Clusters.size(); ++j) {
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
					
				if (DDR < DR) {
					LCLU = j;
					DR = DDR;
				}
			}

			if (DR > RJE) 
				LCLU = -1;
				
			// check if photon  isolated from clusters
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
			Real64_t EDEP = 0.0;
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

			if (EDEP - PT > EDMAX) 
				ISOL = false;
		
			// fill /ISOPHO/ with isolated photon
			if (ISOL) {
				// remove pho-cluster from /CLUSTER/
				if (LCLU >= 0) 
					orecord.Clusters.erase(orecord.Clusters.begin() + LCLU);

				// TODO
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
	// sort KPHO, PPHO
	
	// check with partons
	Int32_t IPHO = 0, IPHISO = 0;
	for (int i=0; i<=NSTOP; ++i) {
		const Particle& part = parts[i];
		
		if (part.type == PT_PHOTON) {
			Real64_t PT = part.pT(); 
			Real64_t ETA = part.getEta();
			Real64_t PHI = part.getPhi();
			Real64_t ENER = 0.0;
			Bool_t ISOL = true;
			
			for (int j=0; j<=NSTOP; ++j) {
				if (abs(parts[j].typeID) <= 21 && i != j && !parts[j].isNeutrino()) 
				{
					Real64_t JPT = parts[j].pT(); 
					Real64_t JETA = parts[j].getEta();
					Real64_t JPHI = parts[j].getPhi();
					
					DDR = sqrt(
						pow(ETA - JETA, 2) + 
						pow(JPHI - PHI, 2)
					);
					
					if (abs(JPHI - PHI) > PI)
						DDR = sqrt(
							pow(ETA - JETA, 2) + 
							pow(abs(JPHI - PHI) - 2*PI, 2)
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
