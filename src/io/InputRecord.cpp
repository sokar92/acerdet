#include "InputRecord.h"
using namespace AcerDet::io;

InputRecord::InputRecord(const vector<Particle>& p) : parts(p) {
}

const vector<Particle>& InputRecord::particles() const { return parts; }

const Particle& InputRecord::getMother(const Particle& p) const {
	return p;
}

const Particle& InputRecord::getDaughter1(const Particle& p) const {
	return p;
}

const Particle& InputRecord::getDaughter2(const Particle& p) const {
	return p;
}
