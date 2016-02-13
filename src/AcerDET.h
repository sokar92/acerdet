#ifndef _AcerDet_Library_Included_
#define _AcerDet_Library_Included_

// ---------------------------------------
// -- Include AcerDET 2.0 Configuration -- 
// ---------------------------------------
#include "conf/Configuration.h"

// ------------------------------------------
// --  Include AcerDET 2.0 core components --
// ------------------------------------------
#include "core/IParticleDataProviderFactory.h"
#include "core/ParticleDataProvider.h"
#include "core/IHistogramManager.h"

// -------------------------------------------
// --  Include AcerDET 2.0 input components --
// -------------------------------------------
#include "io/InputRecord.h"
#include "io/OutputRecord.h"

// ---------------------------------------------
// --  Include AcerDET 2.0 analyse components --
// ---------------------------------------------
#include "analyse/IAnalysePhase.h"

#include <vector>

namespace AcerDet {

	//! Main AcerDET class.
	/**
	 * Program for reconstruction, analyse and storing events statistics from MonteCarlo event generators.
	 */
	class AcerDET {
	private:
		core::IHistogramManager     *histos;            /*!< A global instance of core::IHistogramManager. */
		Bool_t                      histos_initialized; /*!< Indicates whether global instance of core::IHistogramManager is already initialized. */
		core::ParticleDataProvider  partProvider;       /*!< A global instance of core::IParticleDataProvider. */
		std::vector<analyse::IAnalysePhase*> analyse_phases;
		
	public:
	
		/**
		 * Constructor.
		 * Creates a new instance of AcerDET analyzer.
		 * \param conf AcerDET configuration/
		 * \param part instance of global particle data provider.
		 * \param hist instance of global histogram manager.
		 */
		AcerDET(
			const conf::Configuration& conf,
			core::IParticleDataProviderFactory* part,
			core::IHistogramManager* hist );
			
		/**
		 * Destructor.
		 * Performs cleanup.
		 */
		~AcerDET();
		
		/**
		 * Prints out AcerDET basic informations and configuration.
		 */
		void printInfo() const;
		
		/**
		 * Analyse single event.
		 * \param inp input record describing random event.
		 * \param out output record to store analyse results in.
		 * \param weigth event ewigth.
		 */
		void analyseRecord( const io::InputRecord& inp, io::OutputRecord& out, Real64_t weigth );
		
		/**
		 * Prints out event analyse statistics.
		 * Foreach analyse algorithm prints it's statistics.
		 */
		void printResults() const;
		
		/**
		 * Stores all histograms in hard drive.
		 * For ROOT framework requires opened ROOT file before execution.
		 */
		void storeHistograms();
	};

}

#endif
