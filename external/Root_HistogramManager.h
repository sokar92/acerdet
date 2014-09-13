#ifndef _AcerDet_External_Root_HistogramManager_
#define _AcerDet_External_Root_HistogramManager_

#include "../src/core/IHistogramManager.h"
using namespace AcerDet::core;

namespace AcerDet {
	namespace external {
		
		class Root_HistogramManager 
			: public IHistogramManager {
		public:
			void init();
			
			void registerHistogram(
				Int32_t id,
				const string& title,
				Int32_t blocks,
				Real64_t minVal,
				Real64_t maxVal );
			
			void insert( Int32_t id, Real64_t value, Real64_t weigth );
			
			void storeHistograms( const string& file );
		};

	}
}

#endif
