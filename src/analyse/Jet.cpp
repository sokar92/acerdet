#include "Jet.h"
#include <cstdio>

using namespace AcerDet::analyse;

Jet::Jet( const Configuration& config ) :
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
	
	histo_bJets				("Jet: multiplicity", 0.0, 10.0, 10),
	histo_delta_phi			("Jet: delta phi jet-barycentre", -0.5, 0.5, 50),
	histo_delta_eta			("Jet: delta eta jet-barycentre", -0.5, 0.5, 50),
	histo_delta_barycenter	("Jet: delta r jet-barycentre", 0.0, 0.5, 50),
	histo_delta_parton		("Jet: delta r jet-parton", 0.0, 0.5, 50),
	histo_pT_bySum			("Jet: pTjet / SumpTParticle", 0.0, 2.0, 50),
	histo_pT_byPart			("Jet: pTjet / pTparton", 0.0, 2.0, 50)
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
	IEVENT++;
/*	
	// smear clusters energy
	if (KEYSME) {      
		for (int i=0; i<orecord.Clusters.size(); ++i) {
			EECLU = orecord.Clusters[i].pT * cosh(orecord.Clusters[i].eta_rec); // Po co tu liczyc jak zaraz jest nadpisane?
			SIGMA = RESHAD(EECLU, orecord.Clusters[i].eta_rec, CALOTH, orecord.Clusters[i].pT, RCONE);
			
			Real64_t coef = orecord.Clusters[i].pT * (1.0 + SIGMA);
			
			//PXCLU = coef * cos (orecord.Clusters[i].phi_rec);
			//PYCLU = coef * sin (orecord.Clusters[i].phi_rec);
			//PZCLU = coef * sinh(orecord.Clusters[i].eta_rec);
			//EECLU = coef * cosh(orecord.Clusters[i].eta_rec);
	
			Particle pClu;
			pClu.momentum = Vector4(
				cos (orecord.Clusters[i].phi_rec),
				sin (orecord.Clusters[i].phi_rec),
				sinh(orecord.Clusters[i].eta_rec),
				cosh(orecord.Clusters[i].eta_rec)
			);
			
			pClu *= coef;

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
			orecord.Clusters[i].eta_rec = ETACLU;
			orecord.Clusters[i].phi_rec = PHICLU;
			orecord.Clusters[i].pT = PTCLU;
		}
	}
	
	// add nonisolated muons to jets
	for (int i=0; i<NMUOX; ++i) {
		ETA = PMUOX(I,3)
		PHI = PMUOX(I,4)
		PT = PMUOX(I,5)

		DR = 100.0;
		MUCLU = -1;
		for (int j=0; j<orecord.Clusters.size(); ++j) {
			DDR = sqrt(
				pow(ETA - orecord.Clusters[j].eta_rec, 2) +
				pow(PHI - orecord.Clusters[j].phi_rec, 2)
			);
			
			if (abs(PHI - orecord.Clusters[j].phi_rec) > PI)
				DDR = sqrt(
					pow(ETA - orecord.Clusters[j].eta_rec, 2) + 
					pow(abs(PHI - orecord.Clusters[j].phi_rec)-2*PI, 2)
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
				pClu.momentum = muVec * orecord.Clusters[MUCLU].pT + PT * Vec;
				
				//PTCLU = SQRT(PXCLU**2+PYCLU**2)
				//ETACLU = SIGN(LOG((SQRT(PTCLU**2+PZCLU**2)+ABS(PZCLU))/PTCLU),PZCLU) 
				//PHICLU = ANGLE(PXCLU,PYCLU)
				PTCLU = pClu.pT();
				ETACLU = pClu.getEta();
				PHICLU = pClu.getPhi();

				PCLU(MUCLU,3) = ETACLU;
				PCLU(MUCLU,4) = PHICLU;
				PCLU(MUCLU,5) = PTCLU;

				KMUOX(I,5) = 0;
			}
		}
	}
	
	// store accepted jets in \JETALL\ common and flagg in common /CLUSTER/
	for (int i=0; i<orecord.Clusters.size(); ++i) {
		if (orecord.Clusters[i].pT > ETJET && 
			abs(orecord.Clusters[i].eta_rec) < ETAJET) 
		{
			// przenies Cluster -> Jet
			//NJET = NJET + 1
			//DO II from 1 to 5 {          
			//	KJET(NJET,II) = KCLU(I,II) 
			//	PJET(NJET,II) = PCLU(I,II)
			//}

			//KJET(NJET,1) = NJET
			//KCLU(I,5) = 0
		}
	}
	
	// histogram NJET
	CALL HF1(IDENT+1,REAL(NJET),1.0)

	// arrange jets in falling E_T sequence
	
	//DO IDU from 1 to NJET {
	//	ETMAX = 0
	//	DO IMU from 1 to NJET {
	//		IF(KJET(IMU,5) == 0) continue
	//		IF(PJET(IMU,5) < ETMAX) continue
	//		IMAX = IMU
	//		ETMAX = PJET(IMU,5)
	//	}

	//	KJET(IMAX,5) = 0
	//	DO II from 1 to 5 {
	//		KDUM(IDU,II) = KJET(IMAX,II)
	//		PDUM(IDU,II) = PJET(IMAX,II)
	//	}
	//}

	//DO I from 1 to NJET {
	//	DO II from 1 to 5 {
	//		KJET(I,II) = KDUM(I,II)
	//		PJET(I,II) = PDUM(I,II)
	//	}
	//	KJET(I,1) = I
	//	KJET(I,5) = 1
	//}
	
	JetData.sortBy_pT(orecord.Jets);

	// reconstruct baricenter of particles
	for (int i=0; i<orecord.Jets.size(); ++i) {
		Real64_t ETAREC = 0, PTREC = 0, PHIREC = 0;

		for (int j=0; j<parts.size(); ++j) {
			const Particle& part = parts[j];
			 
			if (!part.isStable()) 
				continue;
			
			Real64_t PT = part.pT(), PZ = part.pZ();
			if(PT * PT <= PTLRAT * PZ * PZ) 
				continue; 

			Int32_t KC = part.getKfcomp(); 
			if (KC == 0 || KC == 12 || KC == 14 || KC == 16 || KC == 13 || KC == KFINVS) 
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
			DPHIA = abs(PJET(IJET,2) - PHI); // TODO 

			if (DPHIA > PI) 
				DPHIA = DPHIA - 2*PI;

			if (abs(PJET(IJET,1)) < CALOTH && pow(PJET(IJET,1) - ETA, 2 + pow(DPHIA, 2) > pow(RCONE, 2)) continue;
			if (abs(PJET(IJET,1)) > CALOTH && pow(PJET(IJET,1) - ETA, 2 + pow(DPHIA, 2) > pow(RCONE, 2)) continue;

			PTREC = PTREC + PT;
			ETAREC = ETAREC + ETA * PT;
			PHIREC = PHIREC + PHI * PT;
		}

		ETAREC /= PTREC;
		PHIREC /= PTREC;
		DETR = sqrt(
			pow(ETAREC - PJET(IJET,3), 2) +
			pow(PHIREC - PJET(IJET,4), 2)
		);

		CALL HF1(IDENT+11,ETAREC-PJET(IJET,3),1.0)
		CALL HF1(IDENT+12,PHIREC-PJET(IJET,4),1.0)
		CALL HF1(IDENT+13,DETR,1.0)
		CALL HF1(IDENT+14,PJET(IJET,5)/PTREC,1.0)
	}

	for (int i=0; i<orecord.Jets.size(); ++i) {
		PTREC = 0;
		DETRMIN = RCONE;

		for (int j=6; j<parts.size(); ++j) {
			const Particle& part = parts[j];
			 
			if (part.stateID != 21 || abs(part.typeID) > 10) 
				continue;

			PT = part.pT();
			ETA = part.getEta();
			PHI = part.getPhi();

			DPHIA = abs(PJET(IJET,4) - PHI); 
			if (DPHIA > PI) 
				DPHIA -= 2*PI;

			DETR = sqrt(
				pow(ETA - PJET(IJET,3), 2) +
				pow(DPHIA, 2)
			);
							
			if (DETR < DETRMIN) {
				PTREC = PT;
				DETRMIN = DETR;
			}
		}

		if (PTREC != 0) {
			histo_delta_parton.insert(DETRMIN);
			histo_pT_byPart.insert(PJET(IJET,5) / PTREC);
		}
	}
*/
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
	histo_bJets				.print( true ); //IDENT + 1
	histo_delta_phi			.print( true ); //IDENT + 11
	histo_delta_eta			.print( true ); //IDENT + 12
	histo_delta_barycenter	.print( true ); //IDENT + 13
	histo_delta_parton		.print( true ); //IDENT + 23
	histo_pT_bySum			.print( true ); //IDENT + 14
	histo_pT_byPart			.print( true ); //IDENT + 24
}
