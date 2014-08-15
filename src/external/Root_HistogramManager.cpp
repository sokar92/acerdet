#include "Root_HistogramManager.h"
#include <cstring>
using namespace AcerDet::external;

void Root_HistogramManager::init() {
	HistoManager::getInstance();
}

void Root_HistogramManager::registerHisto(
	const string& name,
	const string& title,
	Int32_t blocks,
	Real64_t minVal,
	Real64_t maxVal ) 
{
	HistoManager::getInstance()
		->addTH1F(name.c_str(), title.c_str(), blocks, minVal, maxVal);
}

void Root_HistogramManager::insert( const string& name, Real64_t value ) {
	HistoManager::getInstance()
		->GetHistoTH1F(name.c_str())
		->Fill(value);
}
