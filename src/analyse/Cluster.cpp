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
	IEVENT++;
	Real64_t PTLRAT = 1.0 / pow(sinh(ETACLU), 2.0);
	
	// Find initiator cell: the one with highest pT of not yet used ones.
	Real64_t ETMAX = 0.0;
	Int32_t ICMAX = -1;
	Real64_t ETA, PHI, THETA, PT;
	
	// compute till maxET > ETINI
	vector<Dum>().swap(orecord.vDum); // wyczysc tmp table
	while (true) {
		
		// find indicator cell with maximum associated energy
		for (int i=0; i<orecord.vCell.size(); ++i) {
			if (orecord.vCell[i].K[4] != 2) 
				continue;

			ICMAX = i;
			ETA = orecord.vCell[i].P[0];
			PHI = orecord.vCell[i].P[1];
			ETMAX = orecord.vCell[i].P[4];
		}
		
		// stop condition - maximum energy is less then required minimum
		if (ETMAX < ETINI)
			break;

		// change state of indicator cell to 'computed'
		orecord.vCell[ICMAX].K[4] = 1;

		// create new cluster from this cell
		Dum newDum;
		newDum.K[2] = 1;						// single hit
		newDum.K[3] = orecord.vCell[i].K[3];	// cell -> cluster
		newDum.K[4] = 1;						// state ok
		
		newDum.P[0] = ETA;
		newDum.P[1] = PHI;
		newDum.P[2] = 0; 
		newDum.P[3] = 0;
		newDum.P[4] = 0;
		
		// Sum up unused cells within required distance of initiator.
		for (int i=0; i<orecord.vCell.size(); ++i) {
			if (orecord.vCell[i].K[4] == 0) 
				continue;
				
			DPHIA = abs(orecord.vCell[i].P[1] - PHI); 
			PHIC = orecord.vCell[i].P[1]; 

			if (DPHIA > PI) 
				PHIC = PHIC + sign(2.0 * PI, PHI)

			if (ABS(ETA) < CALOTH && pow((orecord.vCell[i].P[0]-ETA), 2.0) + pow((PHIC-PHI), 2.0) > pow(RCONE, 2.0)) continue;
			if (ABS(ETA) > CALOTH && pow((orecord.vCell[i].P[0]-ETA), 2.0) + pow((PHIC-PHI), 2.0) > pow(RCONE, 2.0)) continue;
			
			orecord.vCell[i].K[4] = -orecord.vCell[i].K[4]; 
			newDum.K[2]++; // another hit
			newDum.K[3] += orecord.vCell[i].K[3];
			newDum.P[2] += orecord.vCell[i].P[4] * orecord.vCell[i].P[0]; 
			newDum.P[3] += orecord.vCell[i].P[4] * PHIC;
			newDum.P[4] += orecord.vCell[i].P[4]; // sum energy in cluster
			
			i++;
		}
		
		// Reject cluster below minimum ET, else accept. 
		if (newDum.P[4] < ETCLU) {
			// revert changes
			for (int i=0;i<orecord.vCell.size();i++	{
				if (orecord.vCell[i].K[4] < 0) 
					orecord.vCell[i].K[4] = -orecord.vCell[i].K[4]; 
			}
		} else {
			newDum.P[2] /= newDum.P[4]; 
			newDum.P[3] /= newDum.P[4]; 
			
			if (abs(newDum.P[3]) > PI) 
				newDum.P[3] -= sign(2.0 * PI, newDum.P[3]) 

			//mark used cells in orecord.vCell
			for (int i=0;i<orecord.vCell.size();i++	{
				if (orecord.vCell[i].K[4] < 0) 
					orecord.vCell[i].K[4] = 0; 
			}
			
			orecord.vDum.push_back(newDum);
		}
	}
	
	// Arrange clusters in falling ET sequence and store clusters in /CLUSTER/
	for (int i=0; i<orecord.vDum.size(); i++) {
		// orecord.vCluster = sort(dum);
	//	ETMAX = 0
	//	DO IDUM = from 1 to NDUM {						// szukaj maxa
	//		IF(KDUM(IDUM,5) == 0) continue
	//		IF(PDUM(IDUM,5) < ETMAX) continue
	//		IMAX = IDUM
	//		ETMAX = PDUM(IDUM,5)
	//	}
	//	KDUM(IMAX,5) = 0 //sorted state
	//	NCLU = NCLU + 1								// zapisz nowy cluster na outputrecord
	//	DO II = from 1 to 5 {
	//		KCLU(NCLU,II) = KDUM(IMAX,II)
	//		PCLU(NCLU,II) = PDUM(IMAX,II)
	//	}
	//	KCLU(NCLU,5) = 1
	}
	
	// histogram NCLU
	// CALL HF1(IDENT+1,REAL(NCLU),1.0)
	
	// reconstruct baricenter of particles
	for (int ICLU = 0; ICLU < orecord.vCluster.size(); ICLU++) {
		ETAREC = 0;
		PTREC = 0;
		PHIREC = 0;

		const vector<Particle>& parts = irecord.particles();
		for (int i=0;i<parts.size();++i) {
			if (!parts[i].isStable()) 
				continue;
				 
			if (parts[i].pT() <= PTLRAT * parts[i].pZ() * parts[i].pZ()) 
				continue;
 
			KC = kfcomp(parts[i].typeID); 
			if (KC == 0 || KC == 12 || KC == 14 || KC == 16 || KC == 13 || KC == KFINVS) 
				continue;

			Real64_t DETPHI = 0.0;
			Real64_t ETA, PHI, THETA, PT;

			if (KEYFLD && kuchge(part.typeID) != 0) {
				if (part.pT() < PTMIN)
					continue;
					
				Real64_t CHRG = kuchge(part.typeID) / 3.0;
				DETPHI = CHRG * part.foldPhi();
			}
			
			PT = part.pT();
			ETA = part.getEta();
			PHI = saturatePi(part.getPhi() + DETPHI);
			
			// czy nie rownowazne saturatePi?
			DPHIA = abs(orecord.vCluster[ICLU].P[1]-PHI);

			if (DPHIA > PI) 
				DPHIA = DPHIA - 2 * PI;

			if (abs(orecord.vCluster[ICLU].P[0]) < CALOTH && pow((orecord.vCluster[ICLU].P[0]-ETA),2) + pow((DPHIA),2) > pow(RCONE,2)) continue;
			if (abs(orecord.vCluster[ICLU].P[0]) > CALOTH && pow((orecord.vCluster[ICLU].P[0]-ETA),2) + pow((DPHIA),2) > pow(RCONE,2)) continue;

			PTREC += PT;
			ETAREC += ETA * PT;
			PHIREC += PHI * PT;
		}

		ETAREC /= PTREC;
		PHIREC /= PTREC;
		DETR = sqrt( pow((ETAREC-orecord.vCluster[ICLU].P[2]),2) + pow((PHIREC-orecord.vCluste[ICLU].P[3]),2) );
      
		// CALL HF1(IDENT+11, ETAREC-PCLU(ICLU,3), 1.0)
		// CALL HF1(IDENT+12, PHIREC-PCLU(ICLU,4), 1.0)
		// CALL HF1(IDENT+13, DETR, 1.0)
		// CALL HF1(IDENT+14, PCLU(ICLU,5)/PTREC, 1.0)
	}

	for (int ICLU=0; ICLU<orecord.vCluster.size(); ICLU++) {
		PTREC = 0;
		DETRMIN = RCONE;

		for (int i=6; i<parts.size(); i++) { 
			if (parts[i].statusID != 21 || abs(parts[i].typeID) > 10) 
				continue;

			PT = parts[i].pT();
			ETA = parts[i].getEta(); 
			PHI = parts[i].getPhi();
			
			DPHIA = abs(orecord.vCluster[ICLU].P[3] - PHI); 
			if (DPHIA > PI) 
				DPHIA = DPHIA - 2 * PI;
			
			DETR = sqrt( pow((ETA-orecord.vCluster[ICLU].P[2]), 2) + pow(DPHIA, 2))
			IF(DETR < DETRMIN) {
				PTREC = PT;
				DETRMIN = DETR;
			}
		}

		if (PTREC) {
			// CALL HF1(IDENT+23,DETRMIN,1.0)
			// CALL HF1(IDENT+24,PCLU(ICLU,5)/PTREC,1.0)
		}
	}

}
