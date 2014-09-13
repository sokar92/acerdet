#ifndef _AcerDet_Analyse_Cluster_
#define _AcerDet_Analyse_Cluster_

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

		class Cluster {
		private:
			Real64_t ETCLU;
			Real64_t RCONE;
			Real64_t ETACLU;
			Real64_t ETINI;
			
			Real64_t PTMIN;
			Real64_t CALOTH;
			
			Int32_t KEYHID;
			Bool_t KEYFLD;
			Int32_t KFINVS;
			
			Int32_t IEVENT;
			
			IHistogramManager *histoManager;
			Bool_t histoRegistered;
			const ParticleDataProvider& partProvider;

		public:
			Cluster(
				const Configuration&,
				IHistogramManager*,
				const ParticleDataProvider& );

			~Cluster();
			
			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );

			void printResults() const;
		};

	}
}

#endif
