#ifndef _AcerDet_Analyse_Muon_
#define _AcerDet_Analyse_Muon_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/IHistogramManager.h"
#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

namespace AcerDet {
	namespace analyse {

		class Muon {
		private:
			Real64_t ETCLU;
			Real64_t RCONE;
			Real64_t PTMUMIN;
			Real64_t ETAMAX;
			Real64_t RISOLJ;
			Real64_t RDEP;
			Real64_t EDMAX;
			Real64_t PTMUMINT;

			Int32_t KEYHID;
			Bool_t KEYSME;

			Int32_t IEVENT;
			
			//core::Histogram histo_nonisol; //10
			//core::Histogram histo_isol; //11
			//core::Histogram histo_hard; //21
			//core::Histogram histo_sum; //31
		public:
			Muon( const Configuration&, IHistogramManager& );
			~Muon();
			
			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );
			
			void printResults() const;
		};

	}
}

#endif
