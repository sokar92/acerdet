#include "Photon.h"
#include <cstdio>
using namespace AcerDet::analyse;

#include "../core/Smearing.h"
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

void Photon::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {

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

	// variables
	Real64_t PT, ETA, PHI, ENER, JPT, JETA, JPHI, DDR, PTCRU, ENE, EDEP;
		
	// new event to compute
	IEVENT++;

	// reference to particles container
	const vector<Particle>& parts = irecord.particles();

	// look for isolated electrons, sort clusters common
	for (int i=0; i<parts.size(); ++i) {
		const Particle& part = parts[i];
	
		if (part.status != PS_FINAL //|| !part.isStable()
		|| part.pT() == 0) 
			continue;
		
		// analyse photons
		if (part.type == PT_PHOTON) {
			Bool_t ISOL = true;
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
						
				Particle pPho = part;
				pPho.momentum *= (1.0 + Smearing::forPhoton(ENE));
				PT = pPho.pT();
				ENE = pPho.e();
			}

			if (PT < PTLMIN || abs(ETA) > ETAMAX) 
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

				PartData newParton;
				newParton.status = part.statusID;
				newParton.num = i;
				newParton.motherStatus = parts[part.mother].statusID;
				newParton.eta = ETA;
				newParton.phi = PHI;
				newParton.pT = PT;
				
				orecord.Electrons.push_back( newParton );
			}
		}
	}

	// store count in histogram
	histoManager
		->insert(idhist+11, orecord.Photons.size() );

	// arrange photons in falling E_T sequence
	PartData::sortBy_pT( orecord.Photons );
	
	// check with partons
	Int32_t IPHO = 0, IPHOISO = 0;
	for (int i=0; i<parts.size(); ++i) {
		const Particle& part = parts[i];
		// uzyj metody globalnej
// to do
		if (!part.isHardProcess())
			continue;
		
		if (part.type == PT_PHOTON) {
			PT = part.pT(); 
			ETA = part.getEta();
			PHI = part.getPhi();
			ENER = 0.0;
			Bool_t ISOL = true;
// to do			
			for (int j=0; j<parts.size(); ++j) {
				if (parts[j].isHardProcess()
				&& abs(parts[j].pdg_id) <= 21
				&& i != j
				&& !parts[j].isNeutrino()) 
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
