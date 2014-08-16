#include "Root_HistogramManager.h"
#include <cstring>
#include "../../libHistoManager/HistoManager.h"
using namespace AcerDet::external;

void Root_HistogramManager::init() {
	printf ("HISTOS: init()\n");
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
	printf ("HISTOS: registerHistogram(%s, %s)\n", name, title.c_str());
	HistoManager::getInstance()
		->addTH1F(name, title.c_str(), blocks, minVal, maxVal);
}

void Root_HistogramManager::insert( Int32_t id, Real64_t value ) {
	printf ("HISTOS: insert(%d)\n", id);
	HistoManager::getInstance()
		->GetHistoTH1F(id)
		->Fill(value);
}

void Root_HistogramManager::storeHistograms( const string& file ) {
	TFile scanMiniTreeFile(file.c_str(), "recreate");
	scanMiniTreeFile.cd();
	printf ("HISTOS: store(%s)\n", file.c_str());
	HistoManager::getInstance()
		->StoreHistos();
}
