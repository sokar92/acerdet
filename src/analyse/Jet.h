#ifndef _AcerDet_Analyse_Jet_
#define _AcerDet_Analyse_Jet_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/IHistogramManager.h"
#include "../core/ParticleDataProvider.h"
using namespace AcerDet::core;

#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"
using namespace AcerDet::io;

namespace AcerDet {
	namespace analyse {

		class Jet {
		private:
			Real64_t ETJET;
			Real64_t ETAJET;
			Real64_t RCONE;
			Real64_t PTMIN;
			Real64_t CALOTH;

			Int32_t KEYHID;
			Bool_t KEYSME;
			Bool_t KEYFLD;
			Int32_t KFINVS;

			Int32_t IEVENT;
			
			IHistogramManager *histoManager;
			Bool_t histoRegistered;
			const ParticleDataProvider& partProvider;
			
		public:
			Jet(
				const Configuration&,
				IHistogramManager*,
				const ParticleDataProvider& );

			~Jet();
			
			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );
			
			void printResults() const;
		};

	}
}

#endif
