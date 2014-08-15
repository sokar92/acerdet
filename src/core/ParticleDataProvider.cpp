#include "ParticleDataProvider.h"
using namespace AcerDet::core;

template<typename K, typename V>
bool contains(map<K,V> m, const K& key) {
	return m.find(key) != m.end();
}

template<typename K, typename V>
const V& get(map<K,V> m, const K& key) {
	return (*m.find(key)).second;
}

ParticleDataProvider::ParticleDataProvider() {
}

Bool_t ParticleDataProvider::containsInfo(Int32_t pdgId) const {
	return contains(parts, pdgId);
}

void ParticleDataProvider::insertInfo(Int32_t pdgId, Int32_t chgType, Real64_t chg, const string& name) {
	Data data;
	data.chargeType = chgType;
	data.charge = chg;
	data.name = name;
	parts.insert(make_pair(pdgId, data));
}

Int32_t ParticleDataProvider::getChargeType(Int32_t pdgId) const {
	return contains(parts, pdgId) ? get(parts, pdgId).chargeType : 0;
}

Real64_t ParticleDataProvider::getCharge(Int32_t pdgId) const {
	return contains(parts, pdgId) ? get(parts, pdgId).charge : 0.0;
}

const string& ParticleDataProvider::getName(Int32_t pdgId) const {
	string empty = string("");
	return contains(parts, pdgId) ? get(parts, pdgId).name : empty;
}
