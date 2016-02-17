#ifndef _AcerDet_Configuration_
#define _AcerDet_Configuration_

#pragma once

#include <string>
using namespace std;

namespace AcerDet {
	namespace conf {

		//! Global AcerDET configuration.
		/*!
		 * A set of parameters used by analyse algorithms.
		 */
		struct Configuration {
			
			//! Common flags
			/*!
			 * A set of commonly used parameters.
			 */
			struct _Flag {
				int HistogramId;       /*!< id for histograms */
				bool Smearing;         /*!< mearing on = 1, off = 0 */
				bool BField;           /*!< B-field on = 1, off = 0 */
				int SusyParticle;      /*!< code for SUSY LSP particle */

				bool BCJetsLabeling;	/*!< b-jets and c-jets labeling on = 1, off = 0 */
				bool TauJetsLabeling;	/*!< tau-jets labeling on = 1, off = 0 */
				bool JetCalibration;	/*!< jet calibration  on = 1, off = 0 */
				bool Test;		/*!< test histograms  on = 1, off = 0 */
				bool UseFastJet;		/*!< use FastJet clusterization  on = 1, off = 0 */

				/**
				 * Default constructor.
				 * \return new instance of global configuration.
				 */
				_Flag();

				/**
				 * Reads common (globally used) configuration from string.
				 * \param str string to read from.
				 */
				void read( const string& str );
				
				/**
				 * Converts current global configuration into string.
				 * \return converted configuration.
				 */
				string write( ) const;

			} Flag;

			//! Cell-specific flags
			/*!
			 * A set of parameters used by analyse::Cell algorithm.
			 */
			struct _Cell {
				double RapidityCoverage;  /*!< rapidity coverage */
				double MinpT;             /*!< min p_T for B-field */
				double MinEt;             /*!< min E_T for cell */
				double EtaTransition;     /*!< Eta transition in cells granularity.  */
				double GranularityEta;    /*!< Cell eta coordinate granularity. */
				double GranularityPhi;    /*!< Cell phi coordinate granularity. */

				/**
				 * Default constructor.
				 * \return new instance of analyse::Cell configuration.
				 */
				_Cell();

				/**
				 * Reads analyse::Cell configuration from string.
				 * \param str string to read from.
				 */
				void read( const string& str );
				
				/**
				 * Converts current analyse::Cell configuration into string.
				 * \return converted configuration.
				 */
				string write( ) const;

			} Cell;

			//! Cluster-specific flags
			/*!
			 * A set of parameters used by analyse::Cluster algorithm.
			 */
			struct _Cluster {
				double RapidityCoverage;  /*!< rapidity coverage */
				double ConeR;             /*!< cone R for clustering */
				double MinEt;             /*!< minimum E_T for cluster */
				double MinEtInit;         /*!< min E_T for cluster initiator */

				/**
				 * Default constructor.
				 * \return new instance of analyse::Cluster configuration.
				 */
				_Cluster();

				/**
				 * Reads analyse::Cluster configuration from string.
				 * \param str string to read from.
				 */
				void read( const string& str );
				
				/**
				 * Converts current analyse::Cluster configuration into string.
				 * \return converted configuration.
				 */
				string write( ) const;

			} Cluster;

			//! Muon-specific flags
			/*!
			 * A set of parameters used by analyse::Muon algorithm.
			 */
			struct _Muon {
				double MinMomenta;  /*!< minimum muon-momenta to be detected */
				double MaxEta;      /*!< maximum muon eta to be detected */
				double MinIsolRlj;  /*!< min R_lj for muon-isolation */
				double ConeR;       /*!< R_cone for energy deposition */
				double MaxEnergy;   /*!< max energy deposition for isol */

				/**
				 * Default constructor.
				 * \return new instance of analyse::Muon configuration.
				 */
				_Muon();

				/**
				 * Reads analyse::Muon configuration from string.
				 * \param str string to read from.
				 */
				void read( const string& str );
				
				/**
				 * Converts current analyse::Muon configuration into string.
				 * \return converted configuration.
				 */
				string write( ) const;

			} Muon;

			//! Photon-specific flags
			/*!
			 * A set of parameters used by analyse::Photon algorithm.
			 */
			struct _Photon {
				double MinMomenta;  /*!< minimum photon-momenta to be isol */
				double MaxEta;      /*!< maximum photon eta to be isol */
				double MinJetsRlj;  /*!< min R_lj for photon-jet */
				double MinIsolRlj;  /*!< min R_lj for photon-isolation */
				double ConeR;       /*!< R_cone for energy deposition */
				double MaxEnergy;   /*!< max energy deposition for isol */

				/**
				 * Default constructor.
				 * \return new instance of analyse::Photon configuration.
				 */
				_Photon();

				/**
				 * Reads analyse::Electron configuration from string.
				 * \param str string to read from.
				 */
				void read( const string& str );
				
				/**
				 * Converts current analyse::Photon configuration into string.
				 * \return converted configuration.
				 */
				string write( ) const;

			} Photon;

			//! Electron-specific flags
			/*!
			 * A set of parameters used by analyse::Electron algorithm.
			 */
			struct _Electron {
				double MinMomenta;  /*!< minimum electron-momenta to be isol */
				double MaxEta;      /*!< maximum electron eta to be isol */
				double MinJetsRlj;  /*!< min R_lj for electron-jet */
				double MinIsolRlj;  /*!< min R_lj for electron-isolation */
				double ConeR;       /*!< R_cone for energy deposition */
				double MaxEnergy;   /*!< max energy deposition for isol */

				/**
				 * Default constructor.
				 * \return new instance of analyse::Electron configuration.
				 */
				_Electron();

				/**
				 * Reads analyse::Electron configuration from string.
				 * \param str string to read from.
				 */
				void read( const string& str );
				
				/**
				 * Converts current analyse::Electron configuration into string.
				 * \return converted configuration.
				 */
				string write( ) const;

			} Electron;

			//! Jet-specific flags
			/*!
			 * A set of parameters used by analyse::Jet algorithm.
			 */
			struct _Jet {
				double RapidityCoverage;  /*!< rapidity coverage for jets */
				double MinEnergy;         /*!< jets energy_min threshold */

				/**
				 * Default constructor.
				 * \return new instance of analyse::Jet configuration.
				 */
				_Jet();

				/**
				 * Reads analyse::Jet configuration from string.
				 * \param str string to read from.
				 */
				void read( const string& str );
				
				/**
				 * Converts current analyse::Jet configuration into string.
				 * \return converted configuration.
				 */
				string write( ) const;

			} Jet;

			//! BJet-specific flags
			/*!
			 * A set of parameters used by analyse::BJet algorithm.
			 */
			struct _BJet {
				double MinMomenta;  /*!< minimum b-quark pT (after FSR) momenta for b-jet label */
				double MaxEta;      /*!< maximum b-quark eta for b-jet label */
				double MaxRbj;      /*!< max R_bj for b-jet label */

				/**
				 * Default constructor.
				 * \return new instance of analyse::BJet configuration.
				 */
				_BJet();

				/**
				 * Reads analyse::BJet configuration from string.
				 * \param str string to read from.
				 */
				void read( const string& str );
				
				/**
				 * Converts current analyse::BJet configuration into string.
				 * \return converted configuration.
				 */
				string write( ) const;

			} BJet;

			//! CJet-specific flags
			/*!
			 * A set of parameters used by analyse::CJet algorithm.
			 */
			struct _CJet {
				double MinMomenta;  /*!< minimum c-quark pT (after FSR) momenta for c-jet label */
				double MaxEta;      /*!< maximum c-quark eta for c-jet label */
				double MaxRcj;      /*!< max R_cj for c-jet label */

				/**
				 * Default constructor.
				 * \return new instance of analyse::CJet configuration.
				 */
				_CJet();

				/**
				 * Reads analyse::CJet configuration from string.
				 * \param str string to read from.
				 */
				void read( const string& str );
				
				/**
				 * Converts current analyse::CJet configuration into string.
				 * \return converted configuration.
				 */
				string write( ) const;

			} CJet;

			//! Tau-specific flags
			/*!
			 * A set of parameters used by analyse::Tau algorithm.
			 */
			struct _Tau {
				double MinpT;   /*!< minimum tau-had pT for tau-jet label */
				double MaxEta;  /*!< maximum tau-eta for tau-jet label */
				double MinR;    /*!< min R_tauj for tau-jet */
				double MaxR;    /*!< max R_tauj for tau-jet */

				/**
				 * Default constructor.
				 * \return new instance of analyse::Tau configuration.
				 */
				_Tau();

				/**
				 * Reads analyse::Tau configuration from string.
				 * \param str string to read from.
				 */
				void read( const string& str );
				
				/**
				 * Converts current analyse::Tau configuration into string.
				 * \return converted configuration.
				 */
				string write( ) const;

			} Tau;

			//! Mis-specific flags
			/*!
			 * A set of parameters used by analyse::Mis algorithm.
			 */
			struct _Misc {
				double MinEt;  /*!< min E_T for energy in cell to count unused cell */

				/**
				 * Default constructor.
				 * \return new instance of missing energy configuration.
				 */
				_Misc();

				/**
				 * Reads Missing energy configuration from string.
				 * \param str string to read from.
				 */
				void read( const string& str );
				
				/**
				 * Converts current Missing energy configuration into string.
				 * \return converted configuration.
				 */
				string write( ) const;

			} Misc;

			/** Save configuration to file.
			 * \param configuration configuration to save.
			 * \param fileName name of destination file.
			 */
			static void save(
				const Configuration& configuration,
				const std::string& fileName );

			/**
			 * Read configuration from given file.
			 * \param fileName name of file to read from.
			 * \return Configuration read from given file.
			 */
			static Configuration fromFile(
				const std::string& fileName );

			/**
			 * Creates a default Configuration instance. 
			 * \return default configuration
			 */
			static Configuration getDefault();
		};

	}
}

#endif
