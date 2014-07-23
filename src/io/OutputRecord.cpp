#include "OutputRecord.h"
using namespace AcerDet::io;

OutputRecord::OutputRecord() {
}

void OutputRecord::clear() {
	vector<Dum>().swap(vDum);
	vector<Dum>().swap(vCell);

	vector<Particle>().swap(vCluster);
	vector<Particle>().swap(vJet);
	vector<Particle>().swap(vElectron);
	vector<Particle>().swap(vPhoton);
	vector<Particle>().swap(vMuon);
}
