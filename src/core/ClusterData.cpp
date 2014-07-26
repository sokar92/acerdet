#include "ClusterData.h"
using namespace AcerDet::core;

#include <algorithm>
using namespace std;

ClusterData::ClusterData() :
	iepth (0),
	hits (0),
	state (0),
	eta (0.0),
	phi (0.0),
	pT (0.0)
{}

ClusterData::ClusterData(const ClusterData& d) :
	iepth (d.iepth),
	hits (d.hits),
	state (d.state),
	eta (d.eta),
	phi (d.phi),
	pT (d.pT)
{}

ClusterData& ClusterData::operator = (const ClusterData& d) {
	iepth = d.iepth;
	hits = d.hits;
	state = d.state;
	eta = d.eta;
	phi = d.phi;
	pT = d.pT;
	return *this;
}

bool ClusterData::comparator_pT(const ClusterData& l, const ClusterData& r) {
	return l.pT > r.pT;
}

void ClusterData::sortBy_pT(vector<ClusterData>& vec) {
	sort(vec.begin(), vec.end(), comparator_pT);
}
