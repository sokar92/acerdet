#ifndef _AcerDet_Analyse_Electron_
#define _AcerDet_Analyse_Electron_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/IHistogramManager.h"
#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

namespace AcerDet {
	namespace analyse {

		class Electron {
		private:
			Real64_t ETCLU;
			Real64_t RCONE;
			Real64_t PTLMIN;
			Real64_t ETAMAX;
			Real64_t RJE;
			Real64_t RISOLJ;
			Real64_t RDEP;
			Real64_t EDMAX;

			Int32_t KEYHID;
			Bool_t KEYSME;

			Int32_t IEVENT;
			
			//core::Histogram histo_isol;
			//core::Histogram histo_hard;
			//core::Histogram histo_sum;
		public:
			Electron( const Configuration&, IHistogramManager& );
			~Electron();
			
			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );
			
			void printResults() const;
		};

	}
}

#endif
