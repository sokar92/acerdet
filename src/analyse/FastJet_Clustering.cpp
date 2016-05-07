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
	RCONE	( config.Cluster.ConeR ),
	ETACLU	( config.Cluster.RapidityCoverage ),
	ETINI	( config.Cluster.MinEtInit ),
	
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

inline double dAbs(double d) {
	return d < 0 ? -d : d;
}

inline bool closeTo(double a, double b) {
	return dAbs(a - b) < 0.000001;
}

bool cellpXcomparator(const CellData& c1, const CellData& c2) {
	return c1.pX() < c2.pX();
}

void FastJet_Clustering::analyseRecord( const io::InputRecord& irecord, io::OutputRecord& orecord, Real64_t weight ) {
	
	Int32_t idhist = 100 + KEYHID;
	if (!histoRegistered) {
		histoRegistered = true;
		histoManager
			->registerHistogram(idhist+1,  "Cluster: multiplicity", 10, 0.0, 10.0);
	}
	
	IEVENT++;
	
	// FastJet particles container
	vector<PseudoJet> fj_Particles;
	
	for (vector<CellData>::const_iterator it = orecord.Cells.begin(); it != orecord.Cells.end(); it++) {
		const CellData& cell = *it;
		
		PseudoJet jet(cell.pX(), cell.pY(), cell.pZ(), cell.e());
		
		fj_Particles.push_back(jet);
	}
	
	// choose a jet definition
	JetDefinition jet_def(antikt_algorithm, RCONE);
	
	// run the clustering, extract the jets
	ClusterSequence cs(fj_Particles, jet_def);
	
	// fill Cells and Clusters collections in output record
	vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(PTMIN));
	vector<PseudoJet> unusedCells = cs.unclustered_particles();
	vector<PseudoJet> unusedClusters = cs.childless_pseudojets();
	
	// mark all cells in orecord as used
	for (vector<CellData>::iterator it = orecord.Cells.begin(); it != orecord.Cells.end(); it++)
		(*it).status = 0;
	
	// sort cells by pX
	std::sort(orecord.Cells.begin(), orecord.Cells.end(), cellpXcomparator);
	
	// mark unused cells in orecord
	int cells = orecord.Cells.size();
	int changed = 0;
	for (vector<PseudoJet>::const_iterator it = unusedCells.begin(); it != unusedCells.end(); it++) {
		const PseudoJet& jet = *it;
		double jetpX = jet.px();
		
		// find lower bound index
		int l = 0;
		int r = cells - 1;
		while (l != r) {
			int m = (l + r + 1) / 2;
			if (orecord.Cells[m].pX() < jetpX) l = m;
			else r = m - 1;
		}
		int lower_bound = l;
		
		// find upper bound index
		l = 0;
		r = cells - 1;
		while (l != r) {
			int m = (l + r) / 2;
			if (orecord.Cells[m].pX() > jetpX) r = m;
			else l = m + 1;
		}
		int upper_bound = l;
		
		if (lower_bound > upper_bound) {
			printf ("Lower bound is greater than upper bound\n");
			fflush (stdout);
			throw -1;
		}
		
		for (int i = lower_bound; i <= upper_bound; i++) {
			CellData& cell = orecord.Cells[i];
			
			if (cell.status == 0)
				continue;
			
			if (closeTo(cell.pX(), jet.px())
			&& closeTo(cell.pY(), jet.py())
			&& closeTo(cell.pZ(), jet.pz())
			&& closeTo(cell.e(), jet.e()))
			{
				// mark as not used
				cell.status = 2;
				changed++;
				break;
			}
		}
	}
	//printf ("---- changed = %d / %d ------\n", changed, (int)unusedCells.size());
	
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
	
	for (vector<PseudoJet>::const_iterator it = jets.begin(); it != jets.end(); it++) {
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
	
	// fill histogram
	histoManager
	  ->insert(idhist+1, tempClusters.size(), weight);
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
