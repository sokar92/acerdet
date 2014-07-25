#include "ClusterData.h"
using namespace AcerDet::core;

#include <algorithm>
using namespace std;

ClusterData::ClusterData() :
	iepth (0),
	hits (0),
	state (0),
	eta (0.0),
	phi (0.0)
	eT (0.0)
{}

ClusterData::ClusterData(const ClusterData& d) :
	iepth (d.iepth),
	hits (d.hits),
	state (d.state),
	eta (d.eta),
	phi (d.phi),
	eT (d.eT)
{}

ClusterData& ClusterData::operator = (const ClusterData& d) {
	iepth = d.iepth;
	hits = d.hits;
	state = d.state;
	eta = d.eta;
	phi = d.phi;
	eT = d.eT;
	return *this;
}

void ClusterData::sortByE(vector<ClusterData>& vec) {
	return vec;
}
