#include "InputRecord.h"
using namespace AcerDet::io;

/*
 * InputRecord
 */
InputRecord::InputRecord(const vector<Particle>& p) : parts(p) {}

const vector<Particle>& InputRecord::particles() const { return parts; }
