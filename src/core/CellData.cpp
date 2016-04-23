#include "CellData.h"
using namespace AcerDet::core;

#include <algorithm>
#include <cmath>
using namespace std;

CellData::CellData() :
	cellID (-1),
	status (0),
	hits (0),
	eta (0.0),
	phi (0.0),
	pT (0.0)
{}

Real64_t CellData::pX() const { return pT * cos(phi); }

Real64_t CellData::pY() const { return pT * sin(phi); }

Real64_t CellData::pZ() const { return pT * sinh(eta); }

Real64_t CellData::e() const { return pT * cosh(eta); }

bool CellData::comparator_pT(const CellData& l, const CellData& r) {
	return l.pT > r.pT;
}

void CellData::sortBy_pT(vector<CellData>& vec) {
	sort(vec.begin(), vec.end(), comparator_pT);
}
