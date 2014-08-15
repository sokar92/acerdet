#include "Cluster.h"
#include "../core/Functions.h"
#include <cstdio>

using namespace AcerDet::core;
using namespace AcerDet::analyse;

Cluster::Cluster( const Configuration& config, IHistogramManager& histoManager ) :
	ETCLU	( config.Cluster.MinEt ),
	RCONE	( config.Cluster.ConeR ),
	ETACLU	( config.Cluster.RapidityCoverage ),
	ETINI	( config.Cluster.MinEtInit ),
	
	PTMIN	( config.Cell.MinpT ),
	CALOTH	( config.Cell.EtaTransition ),
	
	KEYHID	( config.Flag.HistogramId ),
	KEYFLD	( config.Flag.BField ),
	KFINVS	( config.Flag.SusyParticle ),

	IEVENT	( 0 )//,
	
	//histo_bJets				("Cluster: multiplicity", 0.0, 10.0, 10),
	//histo_delta_phi			("Cluster: delta phi clu-barycentre", -0.5, 0.5, 50),
	//histo_delta_eta			("Cluster: delta eta clu-barycentre", -0.5, 0.5, 50),
	//histo_delta_barycenter	("Cluster: delta r clu-barycentre", 0.0, 0.5, 50),
	//histo_delta_parton		("Cluster: delta r clu-parton", 0.0, 0.5, 50),
	//histo_pT_bySum			("Cluster: pTclu / SumpTParticle", 0.0, 2.0, 50),
	//histo_pT_byPart			("Cluster: pTclu / pTparton", 0.0, 2.0, 50)
{}

Cluster::~Cluster() {}

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
		Real64_t ETMAX = 0.0;
		for (int i=0; i<orecord.Cells.size(); ++i) {
			const CellData& cell = orecord.Cells[i];
			if (cell.state != 2) 
				continue;

			ICMAX = i;
			ETA = cell.eta;
			PHI = cell.phi;
			ETMAX = cell.pT;
		}
		
		// stop condition - maximum energy is less then required minimum
		if (ICMAX < 0 || ETMAX < ETINI) // without first condition Segmentation Fault possible
			break;

		// change state of indicator cell to 'computed'
		orecord.Cells[ICMAX].state = 1; // TODO: create markAsComputed method

		// create new cluster from this cell
		ClusterData newCluster;
		newCluster.cellID = orecord.Cells[ICMAX].cellID; // cell -> cluster
		newCluster.hits = 1;							 // single hit
		newCluster.state = 1;							 // state ok
		
		newCluster.eta = ETA;
		newCluster.phi = PHI;
		newCluster.pT = 0;
		
		// Sum up unused cells within required distance of initiator.
		for (int i=0; i<orecord.Cells.size(); ++i) {
			const CellData& cell = orecord.Cells[i];
			
			if (cell.state == 0) 
				continue;
				
			DPHIA = abs(cell.phi - PHI); 
			PHIC = cell.phi; 

			if (DPHIA > PI) 
				PHIC += sign(2.0 * PI, PHI);

			if (abs(ETA) < CALOTH && pow((cell.eta - ETA), 2.0) + pow((PHIC-PHI), 2.0) > pow(RCONE, 2.0)) continue;
			if (abs(ETA) > CALOTH && pow((cell.eta - ETA), 2.0) + pow((PHIC-PHI), 2.0) > pow(RCONE, 2.0)) continue;
			
			orecord.Cells[i].state = -orecord.Cells[i].state; // negate cell state temporalily 
			newCluster.hits++; // another hit
			newCluster.hits += cell.hits; // ? TODO: makes sense ?
			newCluster.eta_rec += cell.pT * cell.eta; 
			newCluster.phi_rec += cell.pT * PHIC;
			newCluster.pT += cell.pT; // sum energy in cluster
			
			i++;
		}
		
		// Reject cluster below minimum ET, else accept. 
		if (newCluster.pT < ETCLU) {
			// revert changes
			for (int j=0; j<orecord.Cells.size(); j++) {
				if (orecord.Cells[j].state < 0) 
					orecord.Cells[j].state = -orecord.Cells[j].state; 
			}
		} else {
			newCluster.eta_rec /= newCluster.pT; 
			newCluster.phi_rec /= newCluster.pT; 
			
			if (abs(newCluster.phi_rec) > PI) 
				newCluster.phi_rec -= sign(2.0 * PI, newCluster.phi_rec); 

			//mark used cells in orecord.vCell
			for (int j=0; j<orecord.Cells.size(); j++) {
				if (orecord.Cells[j].state < 0) // marked to cluster 
					orecord.Cells[j].state = 0; // joined with cluster
			}
			
			tempClusters.push_back(newCluster);
		}
	}

	// Arrange clusters in falling ET sequence
	ClusterData::sortBy_pT(tempClusters);

	// store in outputrecord
	orecord.Clusters.insert(orecord.Clusters.end(), tempClusters.begin(), tempClusters.end());
	
	// call histogram
	//histo_bJets.insert(tempClusters.size());
	
	// reconstruct baricenter of particles
	const vector<Particle>& parts = irecord.particles();
	for (int ICLU = 0; ICLU < orecord.Clusters.size(); ICLU++) {
		const ClusterData& cluster = orecord.Clusters[ICLU];
		Real64_t ETAREC = 0.0, PTREC = 0.0, PHIREC = 0.0;

		for (int i=0; i<parts.size(); ++i) {
			const Particle& part = parts[i];

			if (!part.isStable()) 
				continue;
			
			Real64_t DETPHI = 0.0;
			Real64_t ETA, PHI, PT, PZ;
			
			PT = part.pT();
			PZ = part.pZ();
			if (PT * PT <= PTLRAT * PZ * PZ) 
				continue;
 
			if (part.type == PT_UNKNOWN || part.isNeutrino() || part.type == PT_MUON || part.type == KFINVS)
				continue;

			/* TODO */
			if (KEYFLD && part.getKuchge() != 0) {
				if (PT < PTMIN)
					continue;
					
				Real64_t CHRG = part.getKuchge() / 3.0;
				DETPHI = CHRG * part.foldPhi();
			}
			
			PT = part.pT();
			ETA = part.getEta();
			PHI = saturatePi(part.getPhi() + DETPHI);
			
			// czy nie rownowazne saturatePi?
			DPHIA = abs(cluster.phi - PHI);

			if (DPHIA > PI) 
				DPHIA -= 2 * PI;

			if (abs(cluster.eta) < CALOTH && pow(cluster.eta - ETA, 2) + pow(DPHIA,2) > pow(RCONE,2)) continue;
			if (abs(cluster.eta) > CALOTH && pow(cluster.eta - ETA, 2) + pow(DPHIA,2) > pow(RCONE,2)) continue;

			PTREC += PT;
			ETAREC += ETA * PT;
			PHIREC += PHI * PT;
		}

		ETAREC /= PTREC;
		PHIREC /= PTREC;
		Real64_t DETR = sqrt( 
			pow((ETAREC - cluster.eta_rec), 2) + 
			pow((PHIREC - cluster.phi_rec), 2) 
		);
		
		// call histograms
		//histo_delta_eta.insert(ETAREC - cluster.eta_rec); //IDENT + 11
		//histo_delta_phi.insert(PHIREC - cluster.phi_rec); //IDENT + 12
		//histo_delta_barycenter.insert(DETR); // IDENT + 13
		//histo_pT_bySum.insert(cluster.pT / PTREC); // IDENT + 14
	}

	for (int ICLU=0; ICLU<orecord.Clusters.size(); ICLU++) {
		const ClusterData& cluster = orecord.Clusters[ICLU];
		Real64_t PTREC = 0.0, DETR;
		Real64_t DETRMIN = RCONE;

		// magic
		for (int i=6; i<parts.size(); i++) { 
			if (parts[i].stateID != 21 || abs(parts[i].typeID) > 10) // TODO: boolean method for this condition
				continue;

			PT = parts[i].pT();
			ETA = parts[i].getEta(); 
			PHI = parts[i].getPhi();
			
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

		if (PTREC) {
			// call histograms
			//histo_delta_parton.insert(DETRMIN); // IDENT + 23
			//histo_pT_byPart.insert(cluster.pT / PTREC); // IDENT + 24
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
	//histo_bJets				.print( true ); //IDENT + 1
	//histo_delta_phi			.print( true ); //IDENT + 11
	//histo_delta_eta			.print( true ); //IDENT + 12
	//histo_delta_barycenter	.print( true ); //IDENT + 13
	//histo_delta_parton		.print( true ); //IDENT + 23
	//histo_pT_bySum			.print( true ); //IDENT + 14
	//histo_pT_byPart			.print( true ); //IDENT + 24
}
