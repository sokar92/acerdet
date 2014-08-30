#ifndef _AcerDet_Analyse_Mis_
#define _AcerDet_Analyse_Mis_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/IHistogramManager.h"
using namespace AcerDet::core;

#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"
using namespace AcerDet::io;

namespace AcerDet {
	namespace analyse {

		class Mis {
		private:
			Real64_t PTMUMIN;
			Real64_t ETAMAX;
			Real64_t ETCELL;
			Real64_t CALOTH;
			
			Int32_t KEYHID;
			Bool_t KEYSME;
			Int32_t KFINVS;
			
			Int32_t IEVENT;
			
			IHistogramManager *histoManager;
			Bool_t histoRegistered;
		
		public:
			Mis( const Configuration&, IHistogramManager* );
			~Mis();
			
			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );
			
			void printResults() const;
		};

	}
}

#endif
