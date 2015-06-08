#include "ClusterData.h"
using namespace AcerDet::core;

#include <algorithm>
using namespace std;

ClusterData::ClusterData() :
	cellID (-1),
	hits (0),
	status (0),
	eta (0.0),
	phi (0.0),
	eta_rec (0.0),
	phi_rec (0.0),
	pT (0.0),
	alreadyUsed (false)
{}

bool ClusterData::comparator_pT(const ClusterData& l, const ClusterData& r) {
	return l.pT > r.pT;
}

void ClusterData::sortBy_pT(vector<ClusterData>& vec) {
	sort(vec.begin(), vec.end(), comparator_pT);
}
