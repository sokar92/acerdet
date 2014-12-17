#include "OutputRecord.h"
using namespace AcerDet::io;

OutputRecord::OutputRecord() {
}

void OutputRecord::clear() {
	vector<CellData>().swap(Cells);
	vector<ClusterData>().swap(Clusters);
	
	vector<PartData>().swap(Muons);
	vector<PartData>().swap(NonisolatedMuons);
	
	vector<PartData>().swap(Electrons);
	vector<PartData>().swap(Photons);
	
	vector<JetData>().swap(Jets);
	Miss.clear();
}
