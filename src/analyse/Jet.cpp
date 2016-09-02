#include "Jet.h"
#include <cstdio>
using namespace AcerDet::analyse;

#include "../core/Smearing.h"
using namespace AcerDet::core;

Jet::Jet(
	const Configuration& config,
	IHistogramManager * histoMng,
	const ParticleDataProvider& partDataProvider ) 
:
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
	histoRegistered( false ),
	partProvider( partDataProvider )	
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

void Jet::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord, Real64_t weight ) {
	
	Int32_t idhist = 500 + KEYHID;
	if (!histoRegistered) {
		histoRegistered = true;
		histoManager
		  ->registerHistogram(idhist+1, "Jet: multiplicity", 30, 0.0, 30);
		histoManager
		  ->registerHistogram(idhist+11, "Jet: delta eta jet-barycentre", 100, -0.5, 0.5);
		histoManager
		  ->registerHistogram(idhist+12, "Jet: delta phi jet-barycentre", 100, -0.5, 0.5);
		histoManager
		  ->registerHistogram(idhist+13, "Jet: delta r jet-barycentre",   100,  0.0, 1.0);
		histoManager
		  ->registerHistogram(idhist+14, "Jet: delta r jet- HP parton",   100,  0.0, 1.0);
		histoManager
		  ->registerHistogram(idhist+23, "Jet: pTjet/pT particles in Rcone",  100,  0.0, 2.0);
		histoManager
		  ->registerHistogram(idhist+24, "Jet: pTjet/pT HP parton",           100,  0.0, 2.0);
	}

	// new event to compute
	IEVENT++;
	
	// reference to particles container
	const vector<Particle>& parts = irecord.particles();
	
	// smear clusters energy
	if (KEYSME) {      
	  for (int i=0; i<orecord.Clusters.size(); ++i) {
	    ClusterData& cluster = orecord.Clusters[i];
	    
	    Real64_t coef = cluster.pT *
	      (1.0 + Smearing::forHadron(cluster.pT * cosh(cluster.eta_rec), cluster.eta_rec, CALOTH));
	    
	    Particle pClu;
	    pClu.momentum = Vector4f(
				     cos (cluster.phi_rec),
				     sin (cluster.phi_rec),
				     sinh(cluster.eta_rec),
				     cosh(cluster.eta_rec)
				     );
	    
	    pClu.momentum *= coef;
	    
	    cluster.eta_rec = -log(max(0.0001, abs(tan(0.5 * pClu.getTheta()))));
	    cluster.phi_rec = pClu.getPhi();
	    cluster.pT      = pClu.pT();
	  }
	}
	
	// add nonisolated muons to jets
	for (int i=0; i<orecord.NonisolatedMuons.size(); ++i) {
	  const ObjectData& muon = orecord.NonisolatedMuons[i];
	  
	  Real64_t DR = 100.0;
	  Int32_t MUCLU = -1;
	  
	  for (int j=0; j<orecord.Clusters.size(); ++j) {
	    const ClusterData& cluster = orecord.Clusters[j];
	    
	    Real64_t DDR = sqrt(
				pow(muon.eta - cluster.eta_rec, 2) +
				pow(muon.phi - cluster.phi_rec, 2)
				);
	    
	    if (abs(muon.phi - cluster.phi_rec) > PI)
	      DDR = sqrt(
			 pow(muon.eta - cluster.eta_rec, 2) + 
			 pow(abs(muon.phi - cluster.phi_rec)-2*PI, 2)
			 );
	    
	    // not chosen yet or less than max
	    if (MUCLU < 0 || DDR < DR) {
	      MUCLU = j;
	      DR = DDR;
	    }
	  }
	  
	  // found candidate
	  if (MUCLU >= 0) {
	    if ((abs(orecord.Clusters[MUCLU].eta_rec) < CALOTH && DR < RCONE) || 
		(abs(orecord.Clusters[MUCLU].eta_rec) > CALOTH && DR < RCONE) 
		) {
	      Vector4f muVec = Vector4f(
					cos (orecord.Clusters[MUCLU].phi_rec),
					sin (orecord.Clusters[MUCLU].phi_rec),
					sinh(orecord.Clusters[MUCLU].eta_rec),
					cosh(orecord.Clusters[MUCLU].eta_rec)
					);
	      Vector4f Vec = Vector4f(cos(muon.phi), sin(muon.phi), sinh(muon.eta), cosh(muon.eta));
	      Particle pClu;
	      pClu.momentum = muVec * orecord.Clusters[MUCLU].pT + Vec * muon.pT;
	      
	      orecord.Clusters[MUCLU].eta_rec = pClu.getEta();
	      orecord.Clusters[MUCLU].phi_rec = pClu.getPhi();
	      orecord.Clusters[MUCLU].pT      = pClu.pT();
	      
	      orecord.NonisolatedMuons[i].alreadyUsed = true;
	    }
	  }
	}
	
	// store accepted jets in \JETALL\ common and flagg in common /CLUSTER/
	for (int i=0; i<orecord.Clusters.size(); ++i) {
	  const ClusterData& cluster = orecord.Clusters[i];
	  
	  // energy > jet_min and angle in range
	  if (cluster.pT > ETJET
	      && abs(cluster.eta_rec) < ETAJET) {
	    JetData newJet;
	    newJet.eta     = cluster.eta;
	    newJet.phi     = cluster.phi;
	    newJet.eta_rec = cluster.eta_rec;
	    newJet.phi_rec = cluster.phi_rec;
	    newJet.pT      = cluster.pT;
	    
	    orecord.Jets.push_back( newJet );
	    orecord.Clusters[i].alreadyUsed = true; // disable converted Cluster
	  }
	}
	
	// histogram NJET
	histoManager
	  ->insert( idhist+1, orecord.Jets.size(), weight );
	
	// arrange jets in falling E_T sequence
	JetData::sortBy_pT( orecord.Jets );
	
	// reconstruct barycenter of particles
	for (int i=0; i<orecord.Jets.size(); ++i) {
	  const JetData& jet = orecord.Jets[i];
	  Real64_t ETAREC = 0, PTREC = 0, PHIREC = 0;
	  
	  for (int j=0; j<parts.size(); ++j) {
	    const Particle& part = parts[j];
	    
	    if (part.status != PS_FINAL) 
	      continue;
	    
	    Real64_t PT = part.pT();
	    Real64_t PZ = part.pZ();
	    Real64_t PTLRAT  = 1.0 / pow(sinh(ETAJET), 2.0);
	    if(PT * PT <= PTLRAT * PZ * PZ) 
	      continue; 
	    
	    if (part.isNeutrino()
		|| part.type == PT_MUON
		|| part.pdg_id == KFINVS)
	      continue;
	    
	    Real64_t DETPHI = 0.0;
	    if (KEYFLD && partProvider.getChargeType(part.pdg_id) != 0) {
	      if (part.pT() < PTMIN) 
		continue;
	      
	      Real64_t CHRG = partProvider.getCharge(part.pdg_id) / 3.0;
	      DETPHI = CHRG * part.foldPhi();
	    }
	    
	    Real64_t PHI = part.getPhi() + DETPHI;
	    Real64_t DPHIA = abs(jet.phi - PHI); 
	    
	    if (DPHIA > PI) 
	      DPHIA = DPHIA - 2*PI;
	    
	    if (abs(jet.eta) < CALOTH
		&& pow(jet.eta - part.getEta(), 2) + pow(DPHIA, 2) > pow(RCONE, 2)) continue;
	    
	    if (abs(jet.eta) > CALOTH
		&& pow(jet.eta - part.getEta(), 2) + pow(DPHIA, 2) > pow(RCONE, 2)) continue;
	    
	    PTREC  += part.pT();
	    ETAREC += part.getEta() * part.pT();
	    PHIREC += PHI * part.pT();
	  }
	  
	  ETAREC /= PTREC;
	  PHIREC /= PTREC;
	  Real64_t DETR = sqrt(
			       pow(ETAREC - jet.eta_rec, 2) +
			       pow(PHIREC - jet.phi_rec, 2)
			       );
	  
	  histoManager
	    ->insert(idhist+11, ETAREC - jet.eta_rec, weight);
	  histoManager
	    ->insert(idhist+12, PHIREC - jet.phi_rec, weight);
	  histoManager
	    ->insert(idhist+13, DETR, weight);
	  histoManager
	    ->insert(idhist+23, jet.pT / PTREC, weight);
	  
	}
	
	// here only hard process particles
	for (int i=0; i<orecord.Jets.size(); ++i) {
	  const JetData& jet = orecord.Jets[i];
	  
	  Real64_t PTREC = 0.0;
	  Real64_t DETRMIN = RCONE;
	  
	  for (int j=6; j<parts.size(); ++j) {
	    const Particle& part = parts[j];
	    
	    if (part.status != PS_HISTORY) 
	      continue;
	    
	    Real64_t DPHIA = abs(part.getPhi() - jet.phi_rec); 
	    if (DPHIA > PI) 
	      DPHIA -= 2*PI;
	    
	    Real64_t DETR = sqrt(
				 pow(part.getEta() - jet.eta_rec, 2) +
				 pow(DPHIA, 2)
				 );
	    
	    if (DETR < DETRMIN) {
	      PTREC = part.pT();
	      DETRMIN = DETR;
	    }
	  }
	  
	  if (PTREC != 0) {
	    histoManager
	      ->insert(idhist+14, DETRMIN, weight);
	    histoManager
	      ->insert(idhist+24, jet.pT / PTREC, weight);
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
}
