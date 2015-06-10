#include "Photon.h"
#include <cstdio>
using namespace AcerDet::analyse;

#include "../core/Smearing.h"
#include "../core/Functions.h"
using namespace AcerDet::core;

Photon::Photon( const Configuration& config, IHistogramManager* histoMng ) :
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
	
	histoManager( histoMng ),
	histoRegistered( false )
	
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

void Photon::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord, Real64_t weigth ) {

	Int32_t idhist = 400 + KEYHID;
	if (!histoRegistered) {
		histoRegistered = true;
		histoManager
			->registerHistogram(idhist+11, "Photon: photon multiplicity ISOLATED", 10, 0.0, 10.0);
		histoManager
			->registerHistogram(idhist+21, "Photon: photon multiplicity HARD", 10, 0.0, 10.0);
		histoManager
			->registerHistogram(idhist+31, "Photon: photon multiplicity HARD+ISOL", 10, 0.0, 10.0);
	}
		
	// new event to compute
	IEVENT++;

	// reference to particles container
	const vector<Particle>& parts = irecord.particles();

	// look for isolated electrons, sort clusters common
	for (int i=0; i<parts.size(); ++i) {
		const Particle& part = parts[i];
	
		if (part.status != PS_FINAL
		|| part.pT() == 0) 
			continue;
		
		// analyse photons
		if (part.type == PT_PHOTON) {

			Bool_t ISOL = true;
			
			// apply smearing
			Real64_t PT = part.pT();
			if (KEYSME) {
				if (parts[i].e() <= 0.0) 
					continue;
						
				Particle pPho = part;
				pPho.momentum *= (1.0 + Smearing::forPhoton(parts[i].e()));
				PT = pPho.pT();
			}			
			

			if (PT < PTLMIN
			|| abs(part.getEta()) > ETAMAX) 
				continue;

			// mark photon-cluster
			Real64_t DR = 100.0;
			Int32_t LCLU = -1;

			for (int j=0; j<orecord.Clusters.size(); ++j) {
				const ClusterData& cluster = orecord.Clusters[j];
				
				Real64_t DDR = sqrt(
					pow(part.getEta() - cluster.eta_rec, 2) +
					pow(part.getPhi() - cluster.phi_rec, 2)
				);
				
				if (abs(part.getPhi() - cluster.phi_rec) > PI)
					DDR = sqrt(
						pow(part.getEta() - cluster.eta_rec, 2) + 
						pow(abs(part.getPhi() - cluster.phi_rec) - 2*PI, 2)
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
				Real64_t DDR = 100.0;
				if (j != LCLU) {
					const ClusterData& cluster = orecord.Clusters[j];
					
					DDR = sqrt(
						pow(part.getEta() - cluster.eta_rec, 2) + 
						pow(part.getPhi() - cluster.phi_rec, 2)
					);
					
					if (abs(part.getPhi() - cluster.phi_rec) > PI)
						DDR = sqrt(
							pow(part.getEta() - cluster.eta_rec, 2) + 
							pow(abs(part.getPhi() - cluster.phi_rec) - 2*PI, 2)
						);
				}

				if (DDR < RISOLJ) 
					ISOL = false;
			}
		
			// check on energy deposition of cells EDEP in cone RDEP
			Real64_t EDEP = 0.0;
			for (int j=0; j<orecord.Cells.size(); ++j) {
				const CellData& cell = orecord.Cells[j];
				
				Real64_t DDR = sqrt(
					pow(part.getEta() - cell.eta, 2) + 
					pow(part.getPhi() - cell.phi, 2)
				);
				
				if (abs(part.getPhi() - cell.phi) > PI)
					DDR = sqrt(
						pow(part.getEta() - cell.eta, 2) + 
						pow(abs(part.getPhi() - cell.phi) - 2*PI, 2)
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

				ObjectData newObject;
				newObject.num = i;
				newObject.pdg_id = part.pdg_id;
				newObject.p = part.momentum;
				newObject.eta = part.getEta();
				newObject.phi = part.getPhi();
				newObject.pT  = PT;
				
				orecord.Photons.push_back( newObject );
			}
		}
	}

	// store count in histogram
	histoManager
		->insert(idhist+11, orecord.Photons.size() );

	// arrange photons in falling E_T sequence
	ObjectData::sortBy_pT( orecord.Photons );
	
	// check with partons
	Int32_t IPHO = 0, IPHOISO = 0;
	for (int i=0; i<parts.size(); ++i) {
		const Particle& part = parts[i];

		if (!isHardProcess(parts, i))
			continue;
		
		if (part.type == PT_PHOTON) { 
			
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

					if (DDR < RISOLJ
					&& parts[j].pT() > ETCLU) 
						ISOL = false;

					if (DDR < RDEP
					&& parts[j].pT() < ETCLU) 
						ENER += parts[j].pT();
				}
			}

			if (ENER > EDMAX)
				ISOL = false;

			if (abs(part.getEta()) < ETAMAX
			&& part.pT() > PTLMIN)
				IPHO++;

			if (abs(part.getEta()) < ETAMAX
			&& part.pT() > PTLMIN && ISOL)
				IPHOISO++;
		}
	}

	// fill histograms
	histoManager
     	->insert(idhist+21, IPHO, 1.0 );
	
	histoManager
     	->insert(idhist+31, IPHOISO, 1.0 );
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
}
