#include "CellData.h"
using namespace AcerDet::core;

#include <algorithm>
using namespace std;

CellData::CellData() :
	iepth (0),
	hits (0),
	state (0),
	eta (0.0),
	phi (0.0)
	pT (0.0)
{}

CellData::CellData(const CellData& d) :
	iepth (d.iepth),
	hits (d.hits),
	state (d.state),
	eta (d.eta),
	phi (d.phi),
	pT (d.pT)
{}

CellData& CellData::operator = (const CellData& d) {
	iepth = d.iepth;
	hits = d.hits;
	state = d.state;
	eta = d.eta;
	phi = d.phi;
	pT = d.pT;
	return *this;
}

void CellData::sortByE(vector<CellData>& vec) {
	return vec;
}
