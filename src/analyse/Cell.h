#ifndef _AcerDet_Analyse_Cell_
#define _AcerDet_Analyse_Cell_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/Histogram.h"
#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

namespace AcerDet {
	namespace analyse {

		class Cell {
		private:
			Real64_t ETACEL;
			Real64_t PTMIN;
			Real64_t ETTHR;
			Real64_t CALOTH;
			Real64_t DBETA;
			Real64_t DBPHI;
			
			Int32_t KEYHID;
			Bool_t KEYFLD;
			Int32_t KFINVS;
			
			Int32_t IEVENT;
			
			core::Histogram histo;
		public:
			Cell( const Configuration& );
			~Cell();
			
			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& ); 
		};

	}
}

#endif
