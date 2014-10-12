#include "CellData.h"
using namespace AcerDet::core;

#include <algorithm>
using namespace std;

CellData::CellData() :
	cellID (-1),
	hits (0),
	status (0),
	eta (0.0),
	phi (0.0),
	pT (0.0)
{}

bool CellData::comparator_pT(const CellData& l, const CellData& r) {
	return l.pT > r.pT;
}

void CellData::sortBy_pT(vector<CellData>& vec) {
	sort(vec.begin(), vec.end(), comparator_pT);
}
