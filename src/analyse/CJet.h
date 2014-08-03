#ifndef _AcerDet_Analyse_C_Jets_
#define _AcerDet_Analyse_C_Jets_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/Histogram.h"
#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

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
			
			core::Histogram histo_cJets; //IDENT + 11
			core::Histogram histo_cQuarks; //IDENT + 21
			core::Histogram histo_delta; //IDENT + 23
			core::Histogram histo_pT; //IDENT + 24
		public:
			CJet( const Configuration& );
			~CJet();

			void printInfo() const;

			void analyseRecord( const io::InputRecord&, io::OutputRecord& );

			void printResults() const;
		};

	}
}

#endif
