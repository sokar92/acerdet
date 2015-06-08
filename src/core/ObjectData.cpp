#include "ObjectData.h"
using namespace AcerDet::core;

#include <algorithm>
using namespace std;

ObjectData::ObjectData() :
	num (0),
	pdg_id (0),
	eta (0.0),
	phi (0.0),
	pT (0.0),
	alreadyUsed (false)
{}

bool ObjectData::comparator_pT(const ObjectData& l, const ObjectData& r) {
	return l.pT > r.pT;
}

void ObjectData::sortBy_pT(vector<ObjectData>& vec) {
	sort(vec.begin(), vec.end(), comparator_pT);
}
