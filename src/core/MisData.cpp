#include "MisData.h"
using namespace AcerDet::core;

MisData::MisData() :
	PXREC (0.0),
	PYREC (0.0),
	PXSUM (0.0),
	PYSUM(0.0),
	PXXCALO (0.0),
	PYYCALO (0.0),
	SUMET(0.0)
{}

void MisData::clear() {
	PXREC = 0.0;
	PYREC = 0.0;
	PXSUM = 0.0;
	PYSUM = 0.0;
	PXXCALO = 0.0;
	PYYCALO = 0.0;
	SUMET = 0.0;
}
