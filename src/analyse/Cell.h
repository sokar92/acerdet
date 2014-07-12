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
			double ETACEL;
			double PTMIN;
			double ETTHR;
			double CALOTH;
			double DBETA;
			double DBPHI;
			
			int KEYHID;
			bool KEYFLD;
			int KFINVS;
			
			int IEVENT;
			
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
