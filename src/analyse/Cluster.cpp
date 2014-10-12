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

void Cluster::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {
	
	Int32_t idhist = 100 + KEYHID;
	if (!histoRegistered) {
		histoRegistered = true;
		histoManager
			->registerHistogram(idhist+1,  "Cluster: multiplicity", 10, 0.0, 10.0);
		histoManager
			->registerHistogram(idhist+11, "Cluster: delta eta clu-barycentre"  , 50, -0.5, 0.5);
		histoManager
			->registerHistogram(idhist+12, "Cluster: delta phi clu-barycentre"  , 50, -0.5, 0.5);
		histoManager
			->registerHistogram(idhist+13, "Cluster: delta r   clu-barycentre"  , 50,  0.0, 0.5);
		histoManager
			->registerHistogram(idhist+23, "Cluster: delta r   clu-parton"      , 50,  0.0, 0.5);
		histoManager
			->registerHistogram(idhist+14, "Cluster: pTclu/SumpTparticle"       , 50,  0.0, 2.0);
		histoManager
			->registerHistogram(idhist+24, "Cluster: pTclu/SumpTparticle"       , 50,  0.0, 2.0);
	}

	// new event to compute
	IEVENT++;
	
	// Find initiator cell: the one with highest pT of not yet used ones.
	Real64_t PTLRAT = 1.0 / pow(sinh(ETACLU), 2.0);
	Real64_t ETA, PHI, PT;
	Real64_t DPHIA, PHIC;

	// temporary clusters container
	vector<ClusterData> tempClusters;
	
	// compute till maxET > ETINI
	while (true) {
		
		// find indicator cell with maximum associated energy
		Int32_t ICMAX = -1;
		for (int i=0; i<orecord.Cells.size(); ++i) {
			const CellData& cell = orecord.Cells[i];
			
			// a cell was used already
			if (cell.status != CellStatus::CREATED) 
				continue;

			if (ICMAX < 0 || cell.pT > orecord.Cells[ICMAX].pT) {
				ICMAX = i;
			}
		}
		
		if (ICMAX < 0  // cell not found (all used) 
		|| orecord.Cells[ICMAX].pT < ETINI) // maximum energy is less then required minimum
			break;

		// change state of indicator cell to 'marked'
		orecord.Cells[ICMAX].status = CellStatus::MARKED;

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
			
			if (cell.status == CellStatus::CLUSTER_JOINED) 
				continue;
				
			DPHIA = abs(cell.phi - orecord.Cells[ICMAX].phi); 
			PHIC = cell.phi; 

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
			
			orecord.Cells[i].status = CellStatus::MARKED;
			newCluster.hits	+= cell.hits; //
			newCluster.eta_rec += cell.pT * cell.eta; // eta_rec = akumulacyjna suma pt * eta 
			newCluster.phi_rec += cell.pT * cell.phi; // phi_rec = akumulacyjna suma pt * phi
			newCluster.pT += cell.pT; // sum energy in cluster
		}
		
		// Reject cluster below minimum ET, else accept. 
		if (newCluster.pT < ETCLU) {
			// revert changes
			for (int j=0; j<orecord.Cells.size(); j++) {
				if (orecord.Cells[j].status == CellStatus::MARKED) 
					orecord.Cells[j].status = CellStatus::CREATED;
			}
			orecord.Cells[ICMAX].status = CellStatus::REJECTED;
		} else {
			newCluster.eta_rec /= newCluster.pT;
			newCluster.phi_rec = saturatePi( newCluster.phi_rec / newCluster.pT );

			//mark used cells in orecord.Cell
			for (int j=0; j<orecord.Cells.size(); j++) {
				if (orecord.Cells[j].status == CellStatus::MARKED) // marked to cluster 
					orecord.Cells[j].status = CellStatus::CLUSTER_JOINED; // joined with cluster
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
		->insert(idhist+1, tempClusters.size());
	
	// reconstruct baricenter of particles
	const vector<Particle>& parts = irecord.particles();
	for (int ICLU = 0; ICLU < orecord.Clusters.size(); ICLU++) {

		const ClusterData& cluster = orecord.Clusters[ICLU];
		Real64_t ETAREC = 0.0, PTREC = 0.0, PHIREC = 0.0;

		for (int i=0; i<parts.size(); ++i) {
			const Particle& part = parts[i];

			//if (!part.isStable())
			if (part.status != PS_FINAL)
				continue;
			
			Real64_t DETPHI = 0.0;
			Real64_t ETA, PHI, PT, PZ;
			
			PT = part.pT();
			PZ = part.pZ();
			
			if (PT * PT <= PTLRAT * PZ * PZ) 
				continue;
 
			if (part.type == PT_UNKNOWN
			|| part.isNeutrino()
			|| part.type == PT_MUON
			|| part.type == KFINVS)
				continue;

			if (KEYFLD && partProvider.getChargeType(part.typeID) != 0) {
				if (PT < PTMIN)
					continue;

				Real64_t CHRG = partProvider.getCharge(part.typeID) / 3.0;
				DETPHI = CHRG * part.foldPhi();
// OBS: charge podaje sensowny bo +-0.33
// printf("DEBUG -> in charge section CHRG = %f DETPHI = %f\n", CHRG, DETPHI);
			}
			
			PT = part.pT();
			ETA = part.getEta();
			PHI = saturatePi(part.getPhi() + DETPHI);
// OBS: PHI oraz cluster.phi grubo poza zakresem! (osiagaja nawet 9) WSZYSTKIE KATY DODATNIE!
// OBS: SPRAWDZ funkcje w Particle (get phi eta ... )
// printf ("DEBUG -> PHI = %f cluster.phi = %f getPhi = %f DETPHI = %f \n", PHI, cluster.phi, part.getPhi(), DETPHI);
			DPHIA = abs(cluster.phi - PHI);

			if (DPHIA > PI) 
				DPHIA -= 2 * PI; // dlaczego nie robic tu saturate ?

			if (abs(cluster.eta) < CALOTH && pow(cluster.eta - ETA, 2) + pow(DPHIA,2) > pow(RCONE,2)) continue;
			if (abs(cluster.eta) > CALOTH && pow(cluster.eta - ETA, 2) + pow(DPHIA,2) > pow(RCONE,2)) continue;

			PTREC += PT;
			ETAREC += ETA * PT;
			PHIREC += PHI * PT;
		}
// OBS: ciezko cokolwiek na temat tego wnioskowac (zmienne akumulacyjne)
// printf ("DEBUG: ETA_rec = %f PHI_rec = %f PTREC = %f \n", ETAREC, PHIREC, PTREC);
		ETAREC /= PTREC;
		PHIREC /= PTREC;
		
		Real64_t DETR = sqrt( 
			pow((ETAREC - cluster.eta_rec), 2) + 
			pow((PHIREC - cluster.phi_rec), 2) 
		);
		
		// fill  histograms
		histoManager
			->insert(idhist + 11, ETAREC - cluster.eta_rec);
		
//		double val = ETAREC - cluster.eta_rec;
//		if (val < -0.5 || 0.5 < val)
//			printf ("cluster_ETA %f - %f = %f out of range [%f, %f]\n", ETAREC, cluster.eta_rec, val, -0.5, 0.5); 
		
		histoManager
			->insert(idhist + 12, PHIREC - cluster.phi_rec); 

// OBS: phirec znacznie poza zakresem [-pi, pi]  osiaga okolo 2pi
//		val = PHIREC - cluster.phi_rec;
//		if (val < -0.5 || 0.5 < val)
//			printf ("cluster_PHI %f - %f = %f out of range [%f, %f]\n", PHIREC, cluster.phi_rec, val, -0.5, 0.5); 	
		
		histoManager
			->insert(idhist + 13, DETR);
// OBS: detr poza zakresem konsekwencja phirec! (patrz definicja DETR)
//		val = DETR;
//		if (val < -0.5 || 0.5 < val)
//			printf ("cluster_DETR %f out of range [%f, %f]\n", val, -0.5, 0.5); 
		
		histoManager
			->insert(idhist + 14, cluster.pT / PTREC);
			
//		val = cluster.pT / PTREC;
//		if (val < 0.0 || 2.0 < val)
//			printf ("cluster_PT %f / %f = %f out of range [%f, %f]\n", cluster.pT, PTREC, val, 0.0, 2.0); 
	}

	for (int ICLU=0; ICLU<orecord.Clusters.size(); ICLU++) {
		const ClusterData& cluster = orecord.Clusters[ICLU];
		Real64_t PTREC = 0.0, DETR;
		Real64_t DETRMIN = RCONE;

		// magic: No, here we want to match to cluster only particles outging
        // magic: from the hard process, like Z->ee, H->gamgam, etc
		// ERW: here only hard-process outgoing particles, so fixed value "6"
        // ERW: should be dropped and start from O BUT this condition should
        // ERW: somehow coded into new flag PT_OutHardProcess
		for (int i=6; i<parts.size(); i++) { 
// OBS: 21 smiedzi!  jaki tu warunek
		//	if (parts[i].statusID != 21 || abs(parts[i].typeID) > 10) // TODO: boolean method for this condition 
		//		continue;
// OBS: nie dochodzi do tej linijki ani razu
printf ("xDebug cluster inside\n");
			PT = parts[i].pT();
			ETA = parts[i].getEta(); 
			PHI = parts[i].getPhi();
	// to samo ? czemu nie saturate		
			DPHIA = abs(cluster.phi_rec - PHI); 
			if (DPHIA > PI) 
				DPHIA = DPHIA - 2 * PI;
			
			DETR = sqrt( 
				pow((ETA - cluster.eta_rec), 2) + 
				pow(DPHIA, 2)
			);
			
			if (DETR < DETRMIN) {
				PTREC = PT;
				DETRMIN = DETR;
			}
		}

		// fill histograms
// OBS: wszystkie ptrec sa rowne 0 dlatego histogramy puste!
		if (PTREC != 0) {
			printf ("xDEBUG %d -> %f\n", IEVENT, PTREC);
			histoManager
				->insert(idhist + 23, DETRMIN);
			
//			if (DETRMIN < 0.0 || 0.5 < DETRMIN)
//				printf ("cluster_DETRMIN %f out of range [%f, %f]\n", DETRMIN, 0.0, 0.5); 
			
			histoManager
				->insert(idhist + 24, cluster.pT / PTREC);
			
//			double val = cluster.pT / PTREC;
//			if (val < 0.0 || 2.0 < val)
//				printf ("cluster_PTREC %f / %f = %f out of range [%f, %f]\n", cluster.pT, PTREC, val, 0.0, 2.0);  
		} //else printf ("DEBUG %d -> %f\n", IEVENT, PTREC);
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
