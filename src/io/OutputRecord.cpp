#include "OutputRecord.h"
using namespace AcerDet::io;

OutputRecord::OutputRecord() {
}

void OutputRecord::clear() {
	vector<CellData>().swap(Cells);
	vector<ClusterData>().swap(Clusters);
	
	vector<PartData>().swap(Muons);
	vector<PartData>().swap(Electrons);
	vector<PartData>().swap(Photons);
}
