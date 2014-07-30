#include "OutputRecord.h"
using namespace AcerDet::io;

OutputRecord::OutputRecord() {
}

void OutputRecord::clear() {
	vector<CellData>().swap(cells);
	vector<ClusterData>().swap(clusters);
}
