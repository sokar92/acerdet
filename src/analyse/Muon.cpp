#include "Muon.h"
#include <cstdio>
using namespace AcerDet::analyse;

#include "../core/Smearing.h"
using namespace AcerDet::core;

Muon::Muon( const Configuration& config, IHistogramManager* histoMng ) :
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
	
	histoManager(histoMng),
	histoRegistered( false )

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
	
	Int32_t idhist = 200 + KEYHID;
	if (!histoRegistered) {
		histoRegistered = true;
		histoManager
			->registerHistogram(idhist+10, "Muon: muon multiplicity NOISOLATED", 10, 0.0, 10.0);
		histoManager
			->registerHistogram(idhist+11, "Muon: muon multiplicity ISOLATED", 10, 0.0, 10.0);
		histoManager
			->registerHistogram(idhist+12, "Muon: muon multiplicity HARD", 10, 0.0, 10.0);
		histoManager
			->registerHistogram(idhist+13, "Muon: muon multiplicity HARD+isol", 10, 0.0, 10.0);
	}

	// new event to compute
	IEVENT++;

	// reference to particles container
	const vector<Particle>& parts = irecord.particles();

	// find last position with '21' status
	Int32_t last21 = -1;
	//Int32_t NSTOP = 0, NSTART = 1;
	for (int i=0; i<parts.size(); ++i) {
		if (parts[i].stateID == 21) {
		//	NSTOP = i-1;
		//	NSTART = i;
			last21 = i;
		}
	}
	
	// look for isolated muons, sort clusters common
	// for (int i=NSTART; i<parts.size(); ++i) {
	for (int i=last21+1; i<parts.size(); ++i) {
		const Particle& part = parts[i];

		if (!part.isStable()) 
			continue;

		if (part.type == PT_MUON) {
			Real64_t ETA, PHI, PT, DDR;

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
				Particle pMuo = part;
				pMuo.momentum /= (1.0 + Smearing::forMuon(PT));
				PT = pMuo.pT();
			}
			
			if (PT < PTMUMIN) 
				continue;
			
			// check if muon is isolated from clusters
			for (int j=0; j<orecord.Clusters.size(); ++j) {
				const ClusterData& cluster = orecord.Clusters[j];
				
				DDR = sqrt( 
					pow(ETA - cluster.eta_rec, 2.0) + 
					pow(PHI - cluster.phi_rec, 2.0) 
				);

				if (abs(PHI - cluster.phi_rec) > PI)
					DDR = sqrt( 
						pow(ETA - cluster.eta_rec, 2.0) + 
						pow(abs(PHI - cluster.phi_rec) - 2.0*PI, 2.0) 
					);

				if (DDR < RISOLJ) 
					ISOL = false;
			}

			// check on energy deposition of cells EDEP in cone RDEP
			Real64_t EDEP = 0.0;
			for (int j=0; j<orecord.Cells.size(); ++j) {
				const CellData& cell = orecord.Cells[j];
				
				DDR = sqrt(
					pow(ETA - cell.eta, 2.0) +
					pow(PHI - cell.phi, 2.0)
				);
				
				if (abs(PHI - cell.phi) > PI)
					DDR = sqrt( 
						pow(ETA - cell.eta, 2) + 
						pow(abs(PHI - cell.phi) - 2*PI, 2) 
					);
				
				if (DDR < RDEP) 
					EDEP += cell.pT;		
			}
			
			if (EDEP > EDMAX) 
				ISOL = false;
			
			PartData newParton;
			newParton.state = part.stateID;
			newParton.particleID = i;
			newParton.motherState = parts[part.mother].stateID;
			newParton.eta = ETA;
			newParton.phi = PHI;
			newParton.pT = PT;
			
			// fill /ISOMUO/ with isolated muon
			if (ISOL) {
				orecord.Muons.push_back(newParton);
			//	NMUO = NMUO + 1
			//	KMUO(NMUO,1) = NMUO				// ktory isol
			//	KMUO(NMUO,2) = K(I,2)			// state z particle
			//	KMUO(NMUO,3) = I				// numer czastki w evencie
			//	KMUO(NMUO,4) = K(K(I,3),2)		// numer matki w evencie
			//	KMUO(NMUO,5) = 1				// stan ?

			//	PMUO(NMUO,1) = ETA
			//	PMUO(NMUO,2) = PHI
			//	PMUO(NMUO,3) = ETA
			//	PMUO(NMUO,4) = PHI
			//	PMUO(NMUO,5) = PT
			
			// fill /NOISOMUO/ with non-isolated muon
			} else {
				orecord.NonisolatedMuons.push_back(newParton);
			//	NMUOX = NMUOX + 1
			//	KMUOX(NMUOX,1) = NMUOX
			//	KMUOX(NMUOX,2) = K(I,2)
          	//	KMUOX(NMUOX,3) = I
			//	KMUOX(NMUOX,4) = K(K(I,3),2)
			//	KMUOX(NMUOX,5) = 1

          	//	PMUOX(NMUOX,1) = ETA
          	//	PMUOX(NMUOX,2) = PHI
          	//	PMUOX(NMUOX,3) = ETA
          	//	PMUOX(NMUOX,4) = PHI
          	//	PMUOX(NMUOX,5) = PT
			}
		}
	}
	
	// call histograms
	histoManager
	  ->insert(idhist + 10,orecord.Muons.size() );
	histoManager
	  ->insert(idhist + 11,orecord.NonisolatedMuons.size() );

	// cross-check with partons
	Int32_t IMUO = 0, IMUOISO = 0;
	//for (int i=0; i<=NSTOP; ++i) {
	for (int i=0; i<=last21; ++i) {
		const Particle& part = parts[i];
		
		if (part.type == PT_MUON) {
			Real64_t PT = part.pT();
			Real64_t ETA = part.getEta(); 
			Real64_t PHI = part.getPhi();
			Real64_t ENER = 0.0;
			Bool_t ISOL = true;
			
			// for (int j=0; j<=NSTOP; ++j) {
			for (int j=0; j<=last21; ++j) {
				if (abs(parts[j].typeID) <= 21
				&& i != j
				&& !parts[j].isNeutrino()) 
				{
					Real64_t JPT = parts[j].pT(); 
					Real64_t JETA = parts[j].getEta(); 
					Real64_t JPHI = parts[j].getPhi();
					
					Real64_t DDR = sqrt(
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
				
			if (abs(ETA) < ETAMAX && PT > PTMUMIN) 
				IMUO++;
				
			if (abs(ETA) < ETAMAX && PT > PTMUMIN && ISOL)
				IMUOISO++;
		}
	}

	// fill histogram
	histoManager
	    ->insert(idhist + 12,  IMUO );
	histoManager
	    ->insert(idhist + 13,  IMUOISO );

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
}
