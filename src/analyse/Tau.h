#ifndef _AcerDet_Analyse_Tau_
#define _AcerDet_Analyse_Tau_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/IHistogramManager.h"
#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

namespace AcerDet {
	namespace analyse {

		class Tau {
		private:
			Real64_t ETJET;
			Real64_t PTTAU;
			Real64_t ETATAU;
			Real64_t RJTAU;
			Real64_t PTFRAC;

			Int32_t KEYHID;
			Bool_t KEYTAU;

			Int32_t IEVENT;
		public:
			Tau( const Configuration&, IHistogramManager& );
			~Tau();

			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );

			void printResults() const;
		};
		
	}
}

#endif
