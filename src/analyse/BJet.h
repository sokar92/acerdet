#ifndef _AcerDet_Analyse_B_Jets_
#define _AcerDet_Analyse_B_Jets_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/Histogram.h"
#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

namespace AcerDet {
	namespace analyse {
		
		class BJet {
		private:
			Real64_t ETJET;
			Real64_t RCONE;
			Real64_t PTBMIN;
			Real64_t ETBMAX;
			Real64_t RJB;

			Int32_t KEYHID;
			Bool_t KEYBCL;

			Int32_t IEVENT;
			
			core::Histogram histo_bJets; //IDENT + 11
			core::Histogram histo_bQuarks; //IDENT + 21
			core::Histogram histo_delta; //IDENT + 23
			core::Histogram histo_pT; //IDENT + 24
		public:
			BJet( const Configuration& );
			~BJet();
			
			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );
			
			void printResults() const;
		};

	}
}

#endif
