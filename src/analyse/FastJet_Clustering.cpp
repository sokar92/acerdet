#include <cstdio>
#include "FastJet_Clustering.h"
#include "../core/Typedefs.h"
#include "../core/Functions.h"
using namespace AcerDet::analyse;

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
using namespace fastjet;

FastJet_Clustering::FastJet_Clustering(
	const Configuration& config,
	IHistogramManager *histoMng,
	const ParticleDataProvider& partDataProvider ) 
:
	ETACEL	( config.Cell.RapidityCoverage ),
	PTMIN	( config.Cell.MinpT ),
	ETTHR	( config.Cell.MinEt ),
	CALOTH	( config.Cell.EtaTransition ),
	DBETA	( config.Cell.GranularityEta ),
	DBPHI	( config.Cell.GranularityPhi ),
	
	ETCLU   ( config.Cluster.MinEt ),
	
	KEYHID	( config.Flag.HistogramId ),
	KEYFLD	( config.Flag.BField ),
	KFINVS	( config.Flag.SusyParticle ),

	IEVENT	( 0 ),
	
	histoManager( histoMng ),
	histoRegistered( false ),
	partProvider( partDataProvider )
{}

FastJet_Clustering::~FastJet_Clustering() {
	histoManager = NULL;
}

void FastJet_Clustering::printInfo() const {
	// print out title
	printf ("***************************************************\n");
	printf ("*                                                 *\n");
	printf ("*     ***************************************     *\n");
	printf ("*     ***   analyse::FastJet_Clustering   ***     *\n");
	printf ("*     ***************************************     *\n");
	printf ("*                                                 *\n");
	printf ("***************************************************\n");

	// print out basic params
	printf ("\n\t clusters definition ...\n");
	printf (" eta coverage  %lf\n", ETACEL);
	printf (" E_T_min cell thresh %lf\n", ETTHR);
	printf (" eta gran. transition %lf\n", CALOTH);
	printf (" gran in eta(central) %lf\n", DBETA);
	printf (" gran in phi(central) %lf\n", DBPHI);
	printf ("\n\t B field apply ....\n");
	printf (" B-field %s\n", KEYFLD ? "on" : "off");
	if (KEYFLD) {
		printf (" p_T min non looping %lf\n", PTMIN);
	}
	printf ("\n");
}

struct ParticleInfo : public PseudoJet::UserInfoBase
{
	Int32_t statusID; // for sure status = PS_FINAL
	ParticleType type;
	Int32_t pdg_id;
			
	Vector4f production;

	Int32_t barcode;
	Int32_t mother;
	pair<Int32_t,Int32_t> daughters;
};

void FastJet_Clustering::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord, Real64_t weight ) {
	
	// reference to particle container
	const vector<Particle>& parts = irecord.particles();
	
	// FastJet particles container
	vector<PseudoJet> fj_Particles;
	
	// Loop over all particles - build input for FastJet
	for (int i=0; i<parts.size(); ++i) {
		const Particle& part = parts[i];

		if (part.status != PS_FINAL)
			continue;
		
		if (part.isNeutrino()
		|| part.type == PT_MUON
		|| part.pdg_id == KFINVS)
			continue;
		
		// set particle momentum
		PseudoJet jet(part.pX(), part.pY(), part.pZ(), part.e());
		
		// set particle additional data
		ParticleInfo info;
		info.statusID = part.statusID;
		info.type = part.type;
		info.pdg_id = part.pdg_id;
		info.production = part.production;
		info.barcode = part.barcode;
		info.mother = part.mother;
		info.daughters = part.daughters;
		
		jet.set_user_info(&info);
		
		// add to collection
		fj_Particles.push_back(jet);
	}
	
	// choose a jet definition
	Real64_t R = 0.7;
	JetDefinition jet_def(antikt_algorithm, R);
	
	// run the clustering, extract the jets
	ClusterSequence cs(fj_Particles, jet_def);
	
	// fill Cells and Clusters collections in output record
	vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(PTMIN));
	vector<PseudoJet> unusedCells = cs.unclustered_particles();
	vector<PseudoJet> unusedClusters = cs.childless_pseudojets();
	
	// fill Cell data
	for (vector<PseudoJet>::const_iterator it = unusedCells.begin(); it != unusedCells.end(); it++) {
		const PseudoJet& jet = *it;
		
		// ignore small cell-impulses
		if (jet.pt() < ETTHR)
			continue;
		
		// fields cellId and hits are ignored - do not need it
		CellData newCell;
		newCell.status = 2;				// CREATED
		newCell.pT = jet.pt();			// pT
		newCell.eta = jet.eta();        // eta
		newCell.phi = jet.phi_std();    // phi \in [-pi, +pi]
		
		// add new cell to output record
		orecord.Cells.push_back(newCell);
	}
	
	// fill Cluster data
	vector<ClusterData> tempClusters;
	for (vector<PseudoJet>::const_iterator it = unusedClusters.begin(); it != unusedClusters.end(); it++) {
		const PseudoJet& jet = *it;
		
		// reject clusters with too small amount of energy
		if (jet.pt() < ETCLU)
			continue;
		
		// fields cellId and hits are ignored - do not need it
		ClusterData newCluster;
		newCluster.status = 1;              // status ok
		newCluster.eta = jet.eta();         // eta
		newCluster.phi = jet.phi_std();     // phi \in [-pi, +pi]
		newCluster.eta_rec = jet.eta();     // eta
		newCluster.phi_rec = jet.phi_std(); // phi \in [-pi, +pi]
		newCluster.pT = jet.pt();
		newCluster.alreadyUsed = false;
		
		tempClusters.push_back(newCluster);
	}
	
	// Arrange clusters in falling ET sequence
	ClusterData::sortBy_pT(tempClusters);

	// store in outputrecord
	orecord.Clusters.insert(orecord.Clusters.end(), tempClusters.begin(), tempClusters.end());
}

void FastJet_Clustering::printResults() const {
	printf ("*************************************************\n");
	printf ("*                                               *\n");
	printf ("*     *************************************     *\n");
	printf ("*     ***         Output from           ***     *\n");
	printf ("*     ***  analyse::FastJet_Clustering  ***     *\n");
	printf ("*     *************************************     *\n");
	printf ("*                                               *\n");
	printf ("*************************************************\n");
	
	printf (" Analysed records: %d\n", IEVENT);
}