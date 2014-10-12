#include "JetData.h"
using namespace AcerDet::core;

#include <algorithm>
using namespace std;

JetData::JetData() :
	type (UNKNOWN),
	eta (0.0),
	phi (0.0),
	eta_rec (0.0),
	phi_rec (0.0),
	pT (0.0)
{}

Bool_t JetData::isBJet() const { return type == B_JET; }
Bool_t JetData::isCJet() const { return type == C_JET; }

bool JetData::comparator_pT(const JetData& l, const JetData& r) {
	return l.pT > r.pT;
}

void JetData::sortBy_pT(vector<JetData>& vec) {
	sort(vec.begin(), vec.end(), comparator_pT);
}
