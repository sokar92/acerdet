#include "Muon.h"
#include <cstdio>
using namespace AcerDet::analyse;

#include "../core/Smearing.h"
#include "../core/Functions.h"
using namespace AcerDet::core;

Muon::Muon( const Configuration& config, IHistogramManager* histoMng ) :
	ETCLU	( config.Cluster.MinEt ),
	RCONE	( config.Cluster.ConeR ),

	PTMUMIN	( config.Muon.MinMomenta ),
	ETAMAX	( config.Muon.MaxEta ),
	RISOLJ	( config.Muon.MinIsolRlj ),
	RDEP	( config.Muon.ConeR ),
	EDMAX	( config.Muon.MaxEnergy ),

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
	printf ("smearing %s\n", KEYSME ? "on" : "off");
	printf ("\n");
}

void Muon::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord, Real64_t weight ) {
	
	Int32_t idhist = 200 + KEYHID;
	if (!histoRegistered) {
		histoRegistered = true;
		histoManager
			->registerHistogram(idhist+10, "Muon: muon multiplicity NOISOLATED", 10, 0.0, 10.0);
		histoManager
			->registerHistogram(idhist+11, "Muon: muon multiplicity ISOLATED", 10, 0.0, 10.0);
		histoManager
			->registerHistogram(idhist+21, "Muon: muon multiplicity HP", 10, 0.0, 10.0);
		histoManager
			->registerHistogram(idhist+31, "Muon: muon multiplicity HP+ISOL", 10, 0.0, 10.0);
	}

	// new event to compute
	IEVENT++;

	// reference to particles container
	const vector<Particle>& parts = irecord.particles();

	// look for isolated muons, sort clusters common
	for (int i=0; i<parts.size(); ++i) {
		const Particle& part = parts[i];

		if (part.status != PS_FINAL)
			continue;

		if (part.type == PT_MUON) {

			// analyse not decayed muons
			Bool_t ISOL = true;
			
			Real64_t PT = part.pT();
			 
			if (abs(part.getEta()) > ETAMAX) 
				continue;  

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
				
				Real64_t DDR = sqrt( 
					pow(part.getEta() - cluster.eta_rec, 2.0) + 
					pow(part.getPhi() - cluster.phi_rec, 2.0) 
				);

				if (abs(part.getPhi() - cluster.phi_rec) > PI)
					DDR = sqrt( 
						pow(part.getEta() - cluster.eta_rec, 2.0) + 
						pow(abs(part.getPhi() - cluster.phi_rec) - 2.0*PI, 2.0) 
					);

				if (DDR < RISOLJ) 
					ISOL = false;
			}

			// check on energy deposition of cells EDEP in cone RDEP
			Real64_t EDEP = 0.0;
			for (int j=0; j<orecord.Cells.size(); ++j) {
				const CellData& cell = orecord.Cells[j];
				
				Real64_t DDR = sqrt(
					pow(part.getEta() - cell.eta, 2.0) +
					pow(part.getPhi() - cell.phi, 2.0)
				);
				
				if (abs(part.getPhi() - cell.phi) > PI)
					DDR = sqrt( 
						pow(part.getEta() - cell.eta, 2) + 
						pow(abs(part.getPhi() - cell.phi) - 2*PI, 2) 
					);
				
				if (DDR < RDEP) 
					EDEP += cell.pT;		
			}
			
			if (EDEP > EDMAX) 
				ISOL = false;

			ObjectData newObject;
			newObject.num = i;
			newObject.pdg_id = part.pdg_id;
			newObject.p = part.momentum;
			newObject.eta = part.getEta();
			newObject.phi = part.getPhi();
			newObject.pT  = PT;
			
			if (ISOL) {
				// fill /ISOMUO/ with isolated muon
				orecord.Muons.push_back(newObject);
			} else {
				// fill /NOISOMUO/ with non-isolated muon
				orecord.NonisolatedMuons.push_back(newObject);
			}
		}
	}
	
	// call histograms
	histoManager
		->insert(idhist + 10, orecord.NonisolatedMuons.size(), weight );
		
	histoManager
		->insert(idhist + 11, orecord.Muons.size(), weight );

	// arrange muons in falling E_T sequence
	ObjectData::sortBy_pT( orecord.Muons );
	ObjectData::sortBy_pT( orecord.NonisolatedMuons );

	// cross-check with partons
	Int32_t IMUO = 0, IMUOISO = 0;
	for (int i=0; i<parts.size(); ++i) {
		const Particle& part = parts[i];
	
		if (!isHardProcess(parts, i))
			continue;
		
		if (part.type == PT_MUON) {
 
			Real64_t ENER = 0.0;
			Bool_t ISOL = true;

			for (int j=0; j<parts.size(); ++j) {
				if (isHardProcess(parts, j)
				&& abs(parts[j].pdg_id) <= 21
				&& i != j
				&& !parts[j].isNeutrino()) 
				{  					
					Real64_t DDR = sqrt(
						pow(part.getEta() - parts[j].getEta(), 2) + 
						pow(part.getPhi() - parts[j].getPhi(), 2) 
					);
					
					if (abs(part.getPhi() - parts[j].getPhi()) > PI)
						DDR = sqrt(
							pow(part.getEta() - parts[j].getEta(), 2) + 
							pow(abs(part.getPhi() - parts[j].getPhi()) - 2*PI, 2) 
						);

					if (DDR < RISOLJ && parts[j].pT() > ETCLU) 
					  ISOL = false;
						
					if (DDR < RDEP && parts[j].pT() < ETCLU) 
					  ENER += parts[j].pT();
				}
			}

			if (ENER > EDMAX) 
			  ISOL = false;
				
			if (abs(part.getEta()) < ETAMAX
                        && part.pT() > PTMUMIN) 
			  IMUO++;
				
			if (abs(part.getEta()) < ETAMAX 
                        && part.pT() > PTMUMIN && ISOL)
			  IMUOISO++;
		}
	}

	// fill histogram
	histoManager
	    ->insert(idhist + 21,  IMUO, weight );
	
	histoManager
	    ->insert(idhist + 31,  IMUOISO, weight );
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
