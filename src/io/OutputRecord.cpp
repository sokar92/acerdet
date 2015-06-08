#include "OutputRecord.h"
using namespace AcerDet::io;

OutputRecord::OutputRecord() {
}

void OutputRecord::clear() {
	vector<CellData>().swap(Cells);
	vector<ClusterData>().swap(Clusters);
	
	vector<ObjectData>().swap(Muons);
	vector<ObjectData>().swap(NonisolatedMuons);
	
	vector<ObjectData>().swap(Electrons);
	vector<ObjectData>().swap(Photons);
	
	vector<JetData>().swap(Jets);
	Miss.clear();
}
