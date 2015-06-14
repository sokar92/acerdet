#include "MissData.h"
using namespace AcerDet::core;

MissData::MissData() :
	PXREC (0.0),
	PYREC (0.0),
	PXSUM (0.0),
	PYSUM (0.0),
	PXXCALO (0.0),
	PYYCALO (0.0),
	PXNUE (0.0),
	PYNUE (0.0),
	SUMET (0.0)
{}

void MissData::clear() {
	PXREC = 0.0;
	PYREC = 0.0;
	PXSUM = 0.0;
	PYSUM = 0.0;
	PXXCALO = 0.0;
	PYYCALO = 0.0;
	PXNUE = 0.0;
	PYNUE = 0.0;
	SUMET = 0.0;
}
