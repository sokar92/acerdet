#ifndef _AcerDet_Analyse_Mis_
#define _AcerDet_Analyse_Mis_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/Histogram.h"
#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

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
			
			core::Histogram histo_reconstructed_pT; //11
			core::Histogram histo_reconstructed_pT_cells; //12
			core::Histogram histo_pTmiss; //13
			core::Histogram histo_pTnu;	//21
		
		public:
			Mis( const Configuration& );
			~Mis();
			
			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );
			
			void printResults() const;
		};

	}
}

#endif
