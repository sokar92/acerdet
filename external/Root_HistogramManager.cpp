#include "Root_HistogramManager.h"
#include <cstring>
#include "libHistoManager/HistoManager.h"
using namespace AcerDet::external;

void Root_HistogramManager::init() {
	HistoManager::getInstance()
		->CreateHistTables();
}

void Root_HistogramManager::registerHistogram(
	Int32_t id,
	const string& title,
	Int32_t blocks,
	Real64_t minVal,
	Real64_t maxVal ) 
{
	char name[256];
	sprintf(name, "hist%.2d", id);
	HistoManager::getInstance()
		->addTH1F(name, title.c_str(), blocks, minVal, maxVal);
}

void Root_HistogramManager::insert( Int32_t id, Real64_t value, Real64_t weigth ) {
	HistoManager::getInstance()
		->GetHistoTH1F(id)
		->Fill(value, weigth);
}

void Root_HistogramManager::storeHistograms() {
	HistoManager::getInstance()
		->StoreHistos();
}
