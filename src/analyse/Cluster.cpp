#include "Cluster.h"
#include "../core/Functions.h"
#include <cstdio>

using namespace AcerDet::core;
using namespace AcerDet::analyse;

Cluster::Cluster( const Configuration& config ) :
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
	
	histo_bJets				("Cluster: multiplicity", 0.0, 10.0, 10),
	histo_delta_phi			("Cluster: delta phi clu-barycentre", -0.5, 0.5, 50),
	histo_delta_eta			("Cluster: delta eta clu-barycentre", -0.5, 0.5, 50),
	histo_delta_barycenter	("Cluster: delta r clu-barycentre", 0.0, 0.5, 50),
	histo_delta_parton		("Cluster: delta r clu-parton", 0.0, 0.5, 50),
	histo_pT_bySum			("Cluster: pTclu / SumpTParticle", 0.0, 2.0, 50),
	histo_pT_byPart			("Cluster: pTclu / pTparton", 0.0, 2.0, 50)
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
		for (int i=0; i<orecord.cells.size(); ++i) {
			if (orecord.cells[i].state != 2) 
				continue;

			ICMAX = i;
			ETA = orecord.cells[i].eta;
			PHI = orecord.cells[i].phi;
			ETMAX = orecord.cells[i].pT;
		}
		
		// stop condition - maximum energy is less then required minimum
		if (ETMAX < ETINI)
			break;

		// change state of indicator cell to 'computed'
		orecord.cells[ICMAX].state = 1; // TODO: create markAsComputed method

		// create new cluster from this cell
		ClusterData newCluster;
		newCluster.iepth = orecord.cells[ICMAX].iepth;	// cell -> cluster
		newCluster.hits = 1;							// single hit
		newCluster.state = 1;							// state ok
		
		newCluster.eta = ETA;
		newCluster.phi = PHI;
		newCluster.pT = 0;
		
		// Sum up unused cells within required distance of initiator.
		for (int i=0; i<orecord.cells.size(); ++i) {
			if (orecord.cells[i].state == 0) 
				continue;
				
			DPHIA = abs(orecord.cells[i].phi - PHI); 
			PHIC = orecord.cells[i].phi; 

			if (DPHIA > PI) 
				PHIC = PHIC + sign(2.0 * PI, PHI);

			if (abs(ETA) < CALOTH && pow((orecord.cells[i].eta - ETA), 2.0) + pow((PHIC-PHI), 2.0) > pow(RCONE, 2.0)) continue;
			if (abs(ETA) > CALOTH && pow((orecord.cells[i].eta - ETA), 2.0) + pow((PHIC-PHI), 2.0) > pow(RCONE, 2.0)) continue;
			
			orecord.cells[i].state = -orecord.cells[i].state; // negate cell state temporalily 
			newCluster.hits++; // another hit
			newCluster.hits += orecord.cells[i].hits;
			// TODO P[2] i P[3] for cluster should be defined
			//newCluster.P[2] += orecord.cells[i].pT * orecord.cells[i].eta; 
			//newCluster.P[3] += orecord.cells[i].pT * PHIC;
			newCluster.pT += orecord.cells[i].pT; // sum energy in cluster
			
			i++;
		}
		
		// Reject cluster below minimum ET, else accept. 
		if (newCluster.pT < ETCLU) {
			// revert changes
			for (int j=0; j<orecord.cells.size(); j++) {
				if (orecord.cells[j].state < 0) 
					orecord.cells[j].state = -orecord.cells[j].state; 
			}
		} else {
			//newCluster.P[2] /= newCluster.pT; 
			//newCluster.P[3] /= newCluster.pT; 
			
			//if (abs(newCluster.P[3]) > PI) 
			//	newCluster.P[3] -= sign(2.0 * PI, newCluster.P[3]) 

			//mark used cells in orecord.vCell
			for (int j=0; j<orecord.cells.size(); j++) {
				if (orecord.cells[j].state < 0) // marked to cluster 
					orecord.cells[j].state = 0; // joined with cluster
			}
			
			tempClusters.push_back(newCluster);
		}
	}

	// Arrange clusters in falling ET sequence
	ClusterData::sortBy_pT(tempClusters);

	// store in /Clusters/
	// TODO : orecord.clusters.insert(tempClusters.begin(), tempClusters.end());
	
	// call histogram
	// .insert(tempClusters.size()) ? ktory to histogram TODO
	
	// reconstruct baricenter of particles
	const vector<Particle>& parts = irecord.particles();
	for (int ICLU = 0; ICLU < orecord.clusters.size(); ICLU++) {
		Real64_t ETAREC = 0.0, PTREC = 0.0, PHIREC = 0.0;

		for (int i=0; i<parts.size(); ++i) {
			const Particle& part = parts[i];

			if (!part.isStable()) 
				continue;
				 
			if (part.pT() <= PTLRAT * part.pZ() * part.pZ()) 
				continue;
 
			Int32_t KC = part.getKfcomp(); 
			if (KC == 0 || KC == 12 || KC == 14 || KC == 16 || KC == 13 || KC == KFINVS) 
				continue;

			Real64_t DETPHI = 0.0;
			Real64_t ETA, PHI, PT;

			if (KEYFLD && part.getKuchge() != 0) {
				if (part.pT() < PTMIN)
					continue;
					
				Real64_t CHRG = part.getKuchge() / 3.0;
				DETPHI = CHRG * part.foldPhi();
			}
			
			PT = part.pT();
			ETA = part.getEta();
			PHI = saturatePi(part.getPhi() + DETPHI);
			
			// czy nie rownowazne saturatePi?
			DPHIA = abs(orecord.clusters[ICLU].phi - PHI);

			if (DPHIA > PI) 
				DPHIA = DPHIA - 2 * PI;

			if (abs(orecord.clusters[ICLU].eta) < CALOTH && pow((orecord.clusters[ICLU].eta - ETA),2) + pow((DPHIA),2) > pow(RCONE,2)) continue;
			if (abs(orecord.clusters[ICLU].eta) > CALOTH && pow((orecord.clusters[ICLU].eta - ETA),2) + pow((DPHIA),2) > pow(RCONE,2)) continue;

			PTREC += PT;
			ETAREC += ETA * PT;
			PHIREC += PHI * PT;
		}

		ETAREC /= PTREC;
		PHIREC /= PTREC;
		//Real64_t DETR = sqrt( pow((ETAREC-orecord.clusters[ICLU].P[2]),2) + pow((PHIREC-orecord.clusters[ICLU].P[3]),2) );
		
		// call histograms ... TODO : ktore?
		// IDENT+11, .insert(ETAREC - orecord.clusters[ICLU].P[2]);
		// IDENT+12, .insert(PHIREC - orecord.clusters[ICLU].P[3]);
		// IDENT+13, .insert(DETR);
		// IDENT+14, .insert(orecord.clusters[ICLU].pT / PTREC);
	}

	for (int ICLU=0; ICLU<orecord.clusters.size(); ICLU++) {
		Real64_t PTREC = 0.0, DETR;
		Real64_t DETRMIN = RCONE;

		// magic
		for (int i=6; i<parts.size(); i++) { 
			if (parts[i].stateID != 21 || abs(parts[i].typeID) > 10) // TODO: boolean method for this condition
				continue;

			PT = parts[i].pT();
			ETA = parts[i].getEta(); 
			PHI = parts[i].getPhi();
			
			//DPHIA = abs(orecord.clusters[ICLU].P[3] - PHI); 
			if (DPHIA > PI) 
				DPHIA = DPHIA - 2 * PI;
			
			//DETR = sqrt( pow((ETA - orecord.clusters[ICLU].P[2]), 2) + pow(DPHIA, 2));
			if (DETR < DETRMIN) {
				PTREC = PT;
				DETRMIN = DETR;
			}
		}

		if (PTREC) {
			// call histograms ... TODO : ktore ?
			// IDENT+23, .insert(DETRMIN);
			// IDENT+24, .insert(orecord.clusters[ICLU].pT / PTREC);
		}
	}
}
