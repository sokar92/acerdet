#ifndef _AcerDet_Analyse_C_Jets_
#define _AcerDet_Analyse_C_Jets_

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

		class CJet {
		private:
			Real64_t ETJET;
			Real64_t RCONE;
			Real64_t PTCMIN;
			Real64_t ETCMAX;
			Real64_t RJC;

			Int32_t KEYHID;
			Bool_t KEYBCL;

			Int32_t IEVENT;

			IHistogramManager *histoManager;
			Bool_t histoRegistered;
			
		public:
			CJet( const Configuration&, IHistogramManager* );
			~CJet();

			void printInfo() const;

			void analyseRecord( const io::InputRecord&, io::OutputRecord& );

			void printResults() const;
		};

	}
}

#endif
