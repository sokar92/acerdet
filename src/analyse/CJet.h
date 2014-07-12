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
			double ETJET;
			double RCONE;
			double PTCMIN;
			double ETCMAX;
			double RJC;

			int KEYHID;
			bool KEYBCL;

			int IEVENT;
			
			core::Histogram histo_bJets; //IDENT + 11
			core::Histogram histo_bQuarks; //IDENT + 21
			core::Histogram histo_delta; //IDENT + 23
			core::Histogram histo_pT; //IDENT + 24
		public:
			CJet( const Configuration& );
			~CJet();
			
			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );
		};

	}
}

#endif
