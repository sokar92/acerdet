#include "PartData.h"
using namespace AcerDet::core;

#include <algorithm>
using namespace std;

PartData::PartData() :
	num (0),
	barcode (-1),
	motherStatus (0),
	status (0),
	eta (0.0),
	phi (0.0),
	pT (0.0)
{}

bool PartData::comparator_pT(const PartData& l, const PartData& r) {
	return l.pT > r.pT;
}

void PartData::sortBy_pT(vector<PartData>& vec) {
	sort(vec.begin(), vec.end(), comparator_pT);
}
