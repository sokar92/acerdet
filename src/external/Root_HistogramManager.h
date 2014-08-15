#ifndef _AcerDet_External_Root_HistogramManager_
#define _AcerDet_External_Root_HistogramManager_

#include "../core/IHistogramManager.h"
using namespace AcerDet::core;

#include "../../libHistoManager/HistoManager.h"

namespace AcerDet {
	namespace external {
		
		class Root_HistogramManager 
			: public IHistogramManager {
		public:
			void init();
			
			void registerHistogram(
				const string& name,
				const string& title,
				Int32_t blocks,
				Real64_t minVal,
				Real64_t maxVal );
			
			void insert( const string& histoName, Real64_t value );
		};

	}
}

#endif
