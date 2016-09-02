#include "Cluster.h"
#include "../core/Functions.h"
#include <cstdio>

using namespace AcerDet::core;
using namespace AcerDet::analyse;

Cluster::Cluster(
	const Configuration& config,
	IHistogramManager *histoMng,
	const ParticleDataProvider& partDataProvider)
:
	ETCLU	( config.Cluster.MinEt ),
	RCONE	( config.Cluster.ConeR ),
	ETACLU	( config.Cluster.RapidityCoverage ),
	ETINI	( config.Cluster.MinEtInit ),
	
	PTMIN	( config.Cell.MinpT ),
	CALOTH	( config.Cell.EtaTransition ),
	
	KEYHID	( config.Flag.HistogramId ),
	KEYFLD	( config.Flag.BField ),
	KFINVS	( config.Flag.SusyParticle ),

	IEVENT	( 0 ),
	
	histoManager(histoMng),
	histoRegistered( false ),
	partProvider( partDataProvider )
{}

Cluster::~Cluster() {
	histoManager = NULL;
}

void Cluster::printInfo() const {
	// print out title
	printf ("****************************************\n");
	printf ("*                                      *\n");
	printf ("*     ****************************     *\n");
	printf ("*     ***   analyse::Cluster   ***     *\n");
	printf ("*     ****************************     *\n");
	printf ("*                                      *\n");
	printf ("****************************************\n");

	// print out basic params
	printf ("\n\t clusters definition ...\n");
	printf (" E_T_min cluster %lf\n", ETCLU);
	printf (" E_T_min cell initia %lf\n", ETINI);
	printf (" R cone %lf\n", RCONE);
	printf (" eta coverage %lf\n", ETACLU);
	printf (" eta gran. transition %lf\n", CALOTH);
	printf ("\n");
}

void Cluster::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord, Real64_t weight ) {
	
	Int32_t idhist = 100 + KEYHID;
	if (!histoRegistered) {
		histoRegistered = true;
		histoManager
			->registerHistogram(idhist+1,  "Cluster: multiplicity", 30, 0.0, 30.0);
		histoManager
			->registerHistogram(idhist+11, "Cluster: delta eta cluster barycentre"  , 100, -0.5, 0.5);
		histoManager
			->registerHistogram(idhist+12, "Cluster: delta phi cluster barycentre"  , 100, -0.5, 0.5);
		histoManager
			->registerHistogram(idhist+13, "Cluster: delta r   cluster barycentre"  , 100,  0.0, 0.5);
		histoManager
			->registerHistogram(idhist+23, "Cluster: delta r   cluster HPparton"    , 100,  0.0, 0.5);
		histoManager
			->registerHistogram(idhist+14, "Cluster: pTclu/SumpT particle"          , 100,  0.0, 2.0);
		histoManager
			->registerHistogram(idhist+24, "Cluster: pTclu/SumpT HP parton"         , 100,  0.0, 2.0);
	}

	// new event to compute
	IEVENT++;
	
	// Find initiator cell: the one with highest pT of not yet used ones.
	Real64_t PTLRAT = 1.0 / pow(sinh(ETACLU), 2.0);

	// temporary clusters container
	vector<ClusterData> tempClusters;
	
	// compute till maxET > ETINI
	while (true) {
		
		// find indicator cell with maximum associated energy
		Int32_t ICMAX = -1;
		for (int i=0; i<orecord.Cells.size(); ++i) {
			const CellData& cell = orecord.Cells[i];
			
			// a cell was used already
			if (cell.status != 2) 
				continue;

			if (ICMAX < 0 || cell.pT > orecord.Cells[ICMAX].pT) {
				ICMAX = i;
			}
		}
		
		if (ICMAX < 0  // cell not found (all used) 
		|| orecord.Cells[ICMAX].pT < ETINI) // maximum energy is less then required minimum
			break;

		// change state of indicator cell to 'marked'
		orecord.Cells[ICMAX].status = 1;

		// create new cluster from this cell
		ClusterData newCluster;
		newCluster.cellID = orecord.Cells[ICMAX].cellID; // cluster inherits cell global ID (related to position in detector)
		newCluster.hits = 1;							 // single hit  tu powinno byc 0 - liczymy w nastepnej petli jeszcze raz
		newCluster.status = 1;							 // state ok
		
		newCluster.eta = orecord.Cells[ICMAX].eta;		// cluster inherits cell eta
		newCluster.phi = orecord.Cells[ICMAX].phi;		// cluster inherits cell phi
		newCluster.pT = 0;
		
		// Sum up unused cells within required distance of initiator.
		for (int i=0; i<orecord.Cells.size(); ++i) {
			const CellData& cell = orecord.Cells[i];
			
			if (cell.status == 0) 
				continue;
				
			Real64_t DPHIA = abs(cell.phi - orecord.Cells[ICMAX].phi); 
			Real64_t PHIC = cell.phi; 

			if (DPHIA > PI) 
				PHIC += sign(2.0 * PI, orecord.Cells[ICMAX].phi);

			if (abs(orecord.Cells[ICMAX].eta) < CALOTH &&
				pow((cell.eta - orecord.Cells[ICMAX].eta), 2.0) +
				pow((cell.phi - orecord.Cells[ICMAX].phi), 2.0)
				> pow(RCONE, 2.0)) continue;
			
			if (abs(orecord.Cells[ICMAX].eta) > CALOTH && 
				pow((cell.eta - orecord.Cells[ICMAX].eta), 2.0) +
				pow((cell.phi - orecord.Cells[ICMAX].phi), 2.0)
				> pow(RCONE, 2.0)) continue;
			
			orecord.Cells[i].status = -orecord.Cells[i].status;
			newCluster.hits	+= cell.hits; //
			newCluster.eta_rec += cell.pT * cell.eta; // eta_rec = akumulacyjna suma pt * eta 
			newCluster.phi_rec += cell.pT * cell.phi; // phi_rec = akumulacyjna suma pt * phi
			newCluster.pT += cell.pT; // sum energy in cluster
		}
		
		// Reject cluster below minimum ET, else accept. 
		if (newCluster.pT < ETCLU) {
			// revert changes
			for (int j=0; j<orecord.Cells.size(); j++) {
				if (orecord.Cells[j].status < 0) 
					orecord.Cells[j].status = -orecord.Cells[j].status;
			}
		} else {
			newCluster.eta_rec /= newCluster.pT;
			newCluster.phi_rec = saturatePi( newCluster.phi_rec / newCluster.pT );

			//mark used cells in orecord.Cell
			for (int j=0; j<orecord.Cells.size(); j++) {
				if (orecord.Cells[j].status < 0) // marked to cluster 
					orecord.Cells[j].status = 0; // joined with cluster
			}
			
			tempClusters.push_back(newCluster);
		}
	}

	// Arrange clusters in falling ET sequence
	ClusterData::sortBy_pT(tempClusters);

	// store in outputrecord
	orecord.Clusters.insert(orecord.Clusters.end(), tempClusters.begin(), tempClusters.end());
	
	// fill histogram
	histoManager
	  ->insert(idhist+1, tempClusters.size(), weight);
	
	// reconstruct baricenter of particles
	const vector<Particle>& parts = irecord.particles();
	for (int ICLU = 0; ICLU < orecord.Clusters.size(); ICLU++) {

		const ClusterData& cluster = orecord.Clusters[ICLU];
		Real64_t ETAREC = 0.0, PTREC = 0.0, PHIREC = 0.0;

		for (int i=0; i<parts.size(); ++i) {
			const Particle& part = parts[i];

			if (part.status != PS_FINAL)
				continue;
			
			Real64_t PT = part.pT();
			Real64_t PZ = part.pZ();
			if (PT * PT <= PTLRAT * PZ * PZ) 
				continue;

			if (part.isNeutrino()
			|| part.type == PT_MUON
			|| part.pdg_id == KFINVS)
				continue;

			Real64_t DETPHI = 0.0;
			if (KEYFLD && partProvider.getChargeType(part.pdg_id) != 0) {
				if (PT < PTMIN)
					continue;

				Real64_t CHRG = partProvider.getCharge(part.pdg_id) / 3.0;
				DETPHI = CHRG * part.foldPhi();
			}
			
			Real64_t PHI = saturatePi(part.getPhi() + DETPHI);
			Real64_t DPHIA = abs(cluster.phi - PHI);
			if (DPHIA > PI) 
				DPHIA -= 2 * PI;

			if (abs(cluster.eta) < CALOTH
			&& pow(cluster.eta - part.getEta(), 2) + pow(DPHIA, 2) > pow(RCONE, 2)) continue;
			
			if (abs(cluster.eta) > CALOTH
			&& pow(cluster.eta - part.getEta(), 2) + pow(DPHIA, 2) > pow(RCONE, 2)) continue;

			PTREC += part.pT();
			ETAREC += part.getEta() * part.pT();
			PHIREC += PHI * part.pT();
		}

		ETAREC /= PTREC;
		PHIREC /= PTREC;
		
		Real64_t DETR = sqrt( 
			pow((ETAREC - cluster.eta_rec), 2) + 
			pow((PHIREC - cluster.phi_rec), 2) 
		);
		
		// fill  histograms
		histoManager
		  ->insert(idhist + 11, ETAREC - cluster.eta_rec, weight);
		
		histoManager
		  ->insert(idhist + 12, PHIREC - cluster.phi_rec, weight); 
		
		histoManager
		  ->insert(idhist + 13, DETR, weight);
		
		histoManager
		  ->insert(idhist + 14, cluster.pT / PTREC, weight);			
	}

	for (int ICLU=0; ICLU<orecord.Clusters.size(); ICLU++) {
		const ClusterData& cluster = orecord.Clusters[ICLU];
		
		Real64_t PTREC = 0.0;
		Real64_t DETRMIN = RCONE;

		// magic: No, here we want to match to cluster only particles outging
        // magic: from the hard process, like Z->ee, H->gamgam, etc
		// ERW: here only hard-process outgoing particles, so fixed value "6"
        // ERW: should be dropped and start from O BUT this condition should
        // ERW: somehow coded into new flag PT_OutHardProcess
		for (int i=6; i<parts.size(); i++) { 
// OBS: 21 smiedzi!  jaki tu warunek
			if (parts[i].statusID != 21 || abs(parts[i].pdg_id) > 10) // TODO: boolean method for this condition 
				continue;
// OBS: nie dochodzi do tej linijki ani razu
			Real64_t DPHIA = abs(cluster.phi_rec - parts[i].getPhi()); 
			if (DPHIA > PI) 
				DPHIA -= 2 * PI;
			
			Real64_t DETR = sqrt( 
				pow(parts[i].getEta() - cluster.eta_rec, 2) + 
				pow(DPHIA, 2)
			);
			
			if (DETR < DETRMIN) {
				PTREC = parts[i].pT();
				DETRMIN = DETR;
			}
		}

		// fill histograms
		if (PTREC != 0) {
			histoManager
			  ->insert(idhist + 23, DETRMIN, weight);
			histoManager
			  ->insert(idhist + 24, cluster.pT / PTREC, weight);
		}
	}
}

void Cluster::printResults() const {
	printf ("**************************************\n");
	printf ("*                                    *\n");
	printf ("*     **************************     *\n");
	printf ("*     ***    Output from     ***     *\n");
	printf ("*     ***  analyse::Cluster  ***     *\n");
	printf ("*     **************************     *\n");
	printf ("*                                    *\n");
	printf ("**************************************\n");
	
	printf (" Analysed records: %d\n", IEVENT);
}
