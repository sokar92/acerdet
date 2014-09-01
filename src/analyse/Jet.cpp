#include "Jet.h"
#include <cstdio>

using namespace AcerDet::analyse;

Jet::Jet( const Configuration& config, IHistogramManager * histoMng ) :
	ETJET	( config.Jet.MinEnergy ),
	ETAJET	( config.Jet.RapidityCoverage ),
	RCONE	( config.Cluster.ConeR ),
	PTMIN	( config.Cell.MinpT ),
	CALOTH	( config.Cell.EtaTransition ),

	KEYHID	( config.Flag.HistogramId ),
	KEYSME	( config.Flag.Smearing ),
	KEYFLD	( config.Flag.BField ),
	KFINVS	( config.Flag.SusyParticle ),

	IEVENT	( 0 ),
	
	histoManager(histoMng),
	histoRegistered( false )
	
{}

Jet::~Jet() {}

void Jet::printInfo() const {
	// print out title
	printf ("************************************\n");
	printf ("*                                  *\n");
	printf ("*     ************************     *\n");
	printf ("*     ***   analyse::Jet   ***     *\n");
	printf ("*     ************************     *\n");
	printf ("*                                  *\n");
	printf ("************************************\n");

	// print out basic params
	printf ("\n\t clusters definition ....\n");
	printf (" R cone %lf\n", RCONE);
	printf ("\t jets definition ....\n");
	printf (" E_T_jets [GeV] %lf\n", ETJET);
	printf (" eta coverage jets %lf\n", ETAJET);
	printf (" smearing %s\n", KEYSME ? "on" : "off");
	printf (" B-field %s\n", KEYFLD ? "on" : "off");
	printf ("\n");
}

void Jet::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord ) {


        int idhist = 500 + KEYHID;

	if (!histoRegistered) {
		histoRegistered = true;
		histoManager
		  ->registerHistogram(idhist+1,"Jet: multiplicity", 50, 0.0, 10);
		histoManager
		  ->registerHistogram(idhist+11,"Jet: delta phi jet-barycentre", 50, -0.5, 0.5);
		histoManager
		  ->registerHistogram(idhist+12,"Jet: delta eta jet-barycentre", 50, -0.5, 0.5);
		histoManager
		  ->registerHistogram(idhist+13,"Jet: delta r jet-barycentre",   50, 0.0, 0.5);
		histoManager
		  ->registerHistogram(idhist+23,"Jet: delta r jet-parton",       50, 0.0, 0.5);
		histoManager
		  ->registerHistogram(idhist+14,"Jet: pTjet/SumpTparticle",      50, 0.0, 2.0);
		histoManager
		  ->registerHistogram(idhist+24,"Jet: pTjet/pTparton",           50, 0.0, 2.0);
	}

	// variables
	Real64_t PT, PZ, ETA, PHI, THETA, DR, DDR, DETR, PTCLU, ETACLU; 
	Real64_t PHICLU, EECLU, SIGMA, DPHIA, DETRMIN, PTREC, PTLRAT, CHRG;

	// new event to compute
	IEVENT++;
	
	// reference to particles container
	const vector<Particle>& parts = irecord.particles();
	
	// smear clusters energy
	if (KEYSME) {      
		for (int i=0; i<orecord.Clusters.size(); ++i) {
			ClusterData& cluster = orecord.Clusters[i];
			
			EECLU = cluster.pT * cosh(cluster.eta_rec); // Po co tu liczyc jak zaraz jest nadpisane?
			//SIGMA = RESHAD(EECLU, cluster.eta_rec, CALOTH, cluster.pT, RCONE);
			
			//Real64_t coef = cluster.pT * (1.0 + SIGMA);
			Real64_t coef = 1.0;
			
			//PXCLU = coef * cos (cluster.phi_rec);
			//PYCLU = coef * sin (cluster.phi_rec);
			//PZCLU = coef * sinh(cluster.eta_rec);
			//EECLU = coef * cosh(cluster.eta_rec);
	
			Particle pClu;
			pClu.momentum = Vector4f(
				cos (cluster.phi_rec),
				sin (cluster.phi_rec),
				sinh(cluster.eta_rec),
				cosh(cluster.eta_rec)
			);
			
			pClu.momentum *= coef;

			THETA = pClu.getTheta();
			//if (abs(PZCLU / EECLU) < 1) {
			//	THETA = ACOS(PZCLU / SQRT(PXCLU**2+PYCLU**2+PZCLU**2))
			//} ELSE {
			//	IF (PZCLU > 0) THETA = 0
			//	IF (PZCLU < 0) THETA = PI
			//}

			ETACLU = -log(max(0.0001, abs(tan(0.5*THETA))));
			//PTCLU  = SQRT(PXCLU**2+PYCLU**2)
			//PHICLU = ANGLE(PXCLU,PYCLU)
			PTCLU = pClu.pT();
			PHICLU = pClu.getPhi();

			//PCLU(I,5) = PTCLU
			//PCLU(I,3) = ETACLU
			//PCLU(I,4) = PHICLU
			cluster.eta_rec = ETACLU;
			cluster.phi_rec = PHICLU;
			cluster.pT = PTCLU;
		}
	}
	
	// add nonisolated muons to jets
	for (int i=0; i<orecord.NonisolatedMuons.size(); ++i) {
		const PartData& muon = orecord.NonisolatedMuons[i];
		
		ETA = muon.eta;
		PHI = muon.phi;
		PT = muon.pT;

		DR = 100.0;
		Int32_t MUCLU = -1;
		for (int j=0; j<orecord.Clusters.size(); ++j) {
			const ClusterData& cluster = orecord.Clusters[j];
			
			DDR = sqrt(
				pow(ETA - cluster.eta_rec, 2) +
				pow(PHI - cluster.phi_rec, 2)
			);
			
			if (abs(PHI - cluster.phi_rec) > PI)
				DDR = sqrt(
					pow(ETA - cluster.eta_rec, 2) + 
					pow(abs(PHI - cluster.phi_rec)-2*PI, 2)
				);

			if (DDR < DR) {
				MUCLU = j;
				DR = DDR;
			}
		}

		// found candidate
		if (MUCLU >= 0) {
			if ((abs(orecord.Clusters[MUCLU].eta_rec) < CALOTH && DR < RCONE) || 
				(abs(orecord.Clusters[MUCLU].eta_rec) > CALOTH && DR < RCONE) 
			) {
				Real64_t coef = orecord.Clusters[MUCLU].pT;
				//PXCLU = PCLU(MUCLU,5) * COS(PCLU(MUCLU,4)) + PT * COS(PHI)
				//PYCLU = PCLU(MUCLU,5) * SIN(PCLU(MUCLU,4)) + PT * SIN(PHI)
				//PZCLU = PCLU(MUCLU,5) * SINH(PCLU(MUCLU,3)) + PT * SINH(ETA)
				//EECLU = PCLU(MUCLU,5) * COSH(PCLU(MUCLU,3)) + PT * COSH(ETA)

				Vector4f muVec = Vector4f(
					cos(orecord.Clusters[MUCLU].phi_rec),
					sin(orecord.Clusters[MUCLU].phi_rec),
					sinh(orecord.Clusters[MUCLU].eta_rec),
					cosh(orecord.Clusters[MUCLU].eta_rec)
				);
				Vector4f Vec = Vector4f(cos(PHI), sin(PHI), sinh(ETA), cosh(ETA));
				Particle pClu;
				pClu.momentum = muVec * orecord.Clusters[MUCLU].pT + Vec * PT;
				
				//PTCLU = SQRT(PXCLU**2+PYCLU**2)
				//ETACLU = SIGN(LOG((SQRT(PTCLU**2+PZCLU**2)+ABS(PZCLU))/PTCLU),PZCLU) 
				//PHICLU = ANGLE(PXCLU,PYCLU)
				PTCLU = pClu.pT();
				ETACLU = pClu.getEta();
				PHICLU = pClu.getPhi();

				orecord.Clusters[MUCLU].eta_rec = ETACLU;
				orecord.Clusters[MUCLU].phi_rec = PHICLU;
				orecord.Clusters[MUCLU].pT = PTCLU;

				// KMUOX(I,5) = 0; ? ze niby co ?
			}
		}
	}
	
	// store accepted jets in \JETALL\ common and flagg in common /CLUSTER/
	for (int i=0; i<orecord.Clusters.size(); ++i) {
		const ClusterData& cluster = orecord.Clusters[i];
		if (cluster.pT > ETJET && abs(cluster.eta_rec) < ETAJET) {
			JetData newJet;
			newJet.eta = cluster.eta;
			newJet.phi = cluster.phi;
			newJet.eta_rec = cluster.eta_rec;
			newJet.phi_rec = cluster.phi_rec;
			newJet.pT = cluster.pT;
			
			orecord.Jets.push_back( newJet );
			// TODO move K

			//KJET(NJET,1) = NJET
			//KCLU(I,5) = 0   // disable converted Cluster
		}
	}
	
	// histogram NJET
	histoManager
  	      ->insert(idhist+1, orecord.Jets.size() );

	// arrange jets in falling E_T sequence
	JetData::sortBy_pT( orecord.Jets );

	// reconstruct baricenter of particles
	for (int i=0; i<orecord.Jets.size(); ++i) {
		const JetData& jet = orecord.Jets[i];
		Real64_t ETAREC = 0, PTREC = 0, PHIREC = 0;

		for (int j=0; j<parts.size(); ++j) {
			const Particle& part = parts[j];
			 
			if (!part.isStable()) 
				continue;
			
			PT = part.pT();
			PZ = part.pZ();
			if(PT * PT <= PTLRAT * PZ * PZ) 
				continue; 

			if (part.type == PT_UNKNOWN || part.isNeutrino() || part.type == PT_MUON || part.type == KFINVS)
				continue;

			Real64_t DETPHI = 0.0;
			if (KEYFLD && part.getKuchge() != 0) {
				PT = part.pT();
				if (PT < PTMIN) 
					continue;

				CHRG = part.getKuchge() / 3.0;
				DETPHI = CHRG * part.foldPhi();
			}

			PT = part.pT();
			ETA = part.getEta();
			PHI = part.getPhi();
			PHI = PHI + DETPHI;
			DPHIA = abs(jet.phi - PHI); 

			if (DPHIA > PI) 
				DPHIA = DPHIA - 2*PI;

			if (abs(jet.eta) < CALOTH && pow(jet.eta - ETA, 2) + pow(DPHIA, 2) > pow(RCONE, 2)) continue;
			if (abs(jet.eta) > CALOTH && pow(jet.eta - ETA, 2) + pow(DPHIA, 2) > pow(RCONE, 2)) continue;

			PTREC = PTREC + PT;
			ETAREC = ETAREC + ETA * PT;
			PHIREC = PHIREC + PHI * PT;
		}

		ETAREC /= PTREC;
		PHIREC /= PTREC;
		DETR = sqrt(
			pow(ETAREC - jet.eta_rec, 2) +
			pow(PHIREC - jet.phi_rec, 2)
		);

		histoManager
		  ->insert(idhist+11, ETAREC - jet.eta_rec);
		histoManager
		  ->insert(idhist+12, PHIREC - jet.phi_rec);
		histoManager
		  ->insert(idhist+13, DETR);
		histoManager
		  ->insert(idhist+23, jet.pT / PTREC);

	}

	for (int i=0; i<orecord.Jets.size(); ++i) {
		const JetData& jet = orecord.Jets[i];
		
		PTREC = 0;
		DETRMIN = RCONE;

		for (int j=6; j<parts.size(); ++j) {
			const Particle& part = parts[j];
			 
			if (part.stateID != 21 || abs(part.typeID) > 10) 
				continue;

			PT = part.pT();
			ETA = part.getEta();
			PHI = part.getPhi();

			DPHIA = abs(jet.phi_rec - PHI); 
			if (DPHIA > PI) 
				DPHIA -= 2*PI;

			DETR = sqrt(
				pow(ETA - jet.eta_rec, 2) +
				pow(DPHIA, 2)
			);
							
			if (DETR < DETRMIN) {
				PTREC = PT;
				DETRMIN = DETR;
			}
		}

		if (PTREC != 0) {
		  histoManager
		     	->insert(idhist+14,DETRMIN);
		  histoManager
			->insert(idhist+24,jet.pT / PTREC);
		}
	}
}

void Jet::printResults() const {
	printf ("**********************************\n");
	printf ("*                                *\n");
	printf ("*     **********************     *\n");
	printf ("*     ***  Output from   ***     *\n");
	printf ("*     ***  analyse::Jet  ***     *\n");
	printf ("*     **********************     *\n");
	printf ("*                                *\n");
	printf ("**********************************\n");
	
	printf (" Analysed records: %d\n", IEVENT);
	//histo_bJets				.print( true ); //IDENT + 1
	//histo_delta_phi			.print( true ); //IDENT + 11
	//histo_delta_eta			.print( true ); //IDENT + 12
	//histo_delta_barycenter	.print( true ); //IDENT + 13
	//histo_delta_parton		.print( true ); //IDENT + 23
	//histo_pT_bySum			.print( true ); //IDENT + 14
	//histo_pT_byPart			.print( true ); //IDENT + 24
}
