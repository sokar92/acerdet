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
/*
	// new event to compute
	IEVENT++;

	// reference to particles container
	const vector<Particle>& parts = irecord.particles();
	Int32_t N = parts.size();

	// znajdz poczatek danych
	Int32_t NSTOP = 0, NSTART = 1;
	for (int i=0; i<parts.size(); ++i) {
		if (parts[i].stateID != 21) {
			NSTOP = i-1;
			NSTART = i;
			break;
		}
	}

	// look for isolated muons, sort clusters common
	for (int i=NSTART; i<parts.size(); ++i) {
		const Particle& part = parts[i];

		if (!parts[i].isStable()) 
			continue;

		if (part.type == PT_MUON) {
			Real64_t ETA, PHI, PT;

			// analyse not decayed muons
			Bool_t ISOL = true;
			
			PT = part.pT();
			if (PT < PTMUMINT) 
				continue;
			
			ETA = part.getEta(); 
			if (abs(ETA) > ETAMAX) 
				continue;  
			
			PHI = part.getPhi();

			// apply smearing
			if (KEYSME) {
				// TODO : SIGMA = RESMUO(PT,ETA,PHI)
				PXMUO = part.pX() / (1.0 + SIGMA);
				PYMUO = part.pY() / (1.0 + SIGMA);
				PZMUO = part.pZ() / (1.0 + SIGMA);
				EEMUO = part.e()  / (1.0 + SIGMA);
				PT    = sqrt(PXMUO*PXMUO + PYMUO*PYMUO);
			}
			
			if (PT < PTMUMIN) 
				continue;
			
			// check if muon is isolated from clusters
			for (int j=0; j<orecord.clusters.size(); ++j) {
				DDR = sqrt( 
					pow(ETA - orecord.clusters[j].eta_rec, 2.0) + 
					pow(PHI - orecord.clusters[j].phi_rec, 2.0) 
				);

				if (abs(PHI - orecord.clusters[j].phi_rec) > PI)
					DDR = sqrt( 
						pow(ETA - orecord.clusters[j].eta_rec, 2.0) + 
						pow(abs(PHI - orecord.clusters[j].phi_rec)-2.0*PI, 2.0) 
					);

				if (DDR < RISOLJ) 
					ISOL = false;
			}

			// check on energy deposition of cells EDEP in cone RDEP
			Real64_t EDEP = 0.0;
			for (int j=0; j<orecord.cells.size(); ++j) {
				DDR = sqrt(
					pow(ETA - orecord.cells[j].eta, 2.0) +
					pow(PHI - orecord.cells[j].phi, 2.0)
				);
				
				if (abs(PHI - orecord.cells[j].phi) > PI)
					DDR = sqrt( 
						pow(ETA - orecord.cells[j].eta, 2) + 
						pow(abs(PHI - orecord.cells[j].phi) - 2*PI, 2) 
					);
				
				if (DDR < RDEP) 
					EDEP += orecord.cells[j].pT;		
			}
			
			if (EDEP > EDMAX) 
				ISOL = false;
// TODO : PARTICLE STRUCTURE				
			// fill /ISOMUO/ with isolated muon
			if (ISOL) {
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
			} else {
			// fill /NOISOMUO/ with non-isolated muon
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
			}
		}
	}
	
	// call histograms
	// IDENT+11, .insert(NMUO);
	// IDENT+10, .insert(NMUOX);
	
	// cross-check with partons
	Int32_t IMUO = 0, IMUOISO = 0;
	for (int i=0; i<=NSTOP; ++i) {
		if (part.type == PT_MUON) {
			PT = part.pT();
			ETA = part.getEta(); 
			PHI = part.getPhi();
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
					
					DDR = sqrt( pow(ETA-JETA, 2) + pow(JPHI-PHI, 2) );
					if (abs(JPHI-PHI) > PI)
						DDR = sqrt( pow(ETA-JETA, 2) + pow(abs(JPHI-PHI)-2*PI, 2) );

					if (DDR < RISOLJ && JPT > ETCLU) 
						ISOL = false;
						
					if (DDR < RDEP && JPT < ETCLU) 
						ENER += JPT;
				}
			}

			if (ENER > EDMAX) 
				ISOL = false;
				
			if (abs(ETA) < ETAMAX && PT > PTMUMIN) 
				IMUO++;
				
			if (abs(ETA) < ETAMAX && PT > PTMUMIN && ISOL)
				IMUOISO++;
		}
	}*/

	// IDENT+21 .insert(IMUO);
	// IDENT+31 .insert(IMUOISO);
}

void Muon::printResults() const {
	printf ("***********************************\n");
	printf ("*                                 *\n");
	printf ("*     ***********************     *\n");
	printf ("*     ***   Output from   ***     *\n");
	printf ("*     ***  analyse::Muon  ***     *\n");
	printf ("*     ***********************     *\n");
	printf ("*                                 *\n");
	printf ("***********************************\n");
	
	printf (" Analysed records: %d\n", IEVENT);
	histo_nonisol	.print( true );
	histo_isol		.print( true );
	histo_hard		.print( true );
	histo_sum		.print( true );
}
