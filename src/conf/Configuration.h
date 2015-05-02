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
				int HistogramId;       /*!< detailed description  */
				bool Smearing;         /*!< Should AcerDET use energy smearing? */
				bool BField;           /*!< detailed description  */
				int SusyParticle;      /*!< PDGid of SUSY particle. */

				bool BCJetsLabeling;	/*!< Should AcerDET label b-jets and c-jets? */
				bool TauJetsLabeling;	/*!< Should AcerDET label tau jets? */
				bool JetCalibration;	/*!< Should AcerDET use calibration for jets? */

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
				double RapidityCoverage;  /*!< detailed description  */
				double MinpT;             /*!< Minimal particle energy to treat it as interesting object. */
				double MinEt;             /*!< Minimal energy to mark parton as a Cell. */
				double EtaTransition;     /*!< detailed description  */
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
				double RapidityCoverage;  /*!< detailed description  */
				double ConeR;             /*!< detailed description  */
				double MinEt;             /*!< detailed description  */
				double MinEtInit;         /*!< detailed description  */

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
				double MinMomenta;  /*!< detailed description  */
				double MaxEta;      /*!< detailed description  */
				double MinIsolRlj;  /*!< detailed description  */
				double ConeR;       /*!< detailed description  */
				double MaxEnergy;   /*!< detailed description  */

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
				double MinMomenta;  /*!< detailed description  */
				double MaxEta;      /*!< detailed description  */
				double MinJetsRlj;  /*!< detailed description  */
				double MinIsolRlj;  /*!< detailed description  */
				double ConeR;       /*!< detailed description  */
				double MaxEnergy;   /*!< detailed description  */

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
				double MinMomenta;  /*!< detailed description  */
				double MaxEta;      /*!< detailed description  */
				double MinJetsRlj;  /*!< detailed description  */
				double MinIsolRlj;  /*!< detailed description  */
				double ConeR;       /*!< detailed description  */
				double MaxEnergy;   /*!< detailed description  */

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
				double RapidityCoverage;  /*!< detailed description  */
				double MinEnergy;         /*!< detailed description  */

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
				double MinMomenta;  /*!< detailed description  */
				double MaxEta;      /*!< detailed description  */
				double MaxRbj;      /*!< detailed description  */

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
				double MinMomenta;  /*!< detailed description  */
				double MaxEta;      /*!< detailed description  */
				double MaxRcj;      /*!< detailed description  */

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
				double MinpT;   /*!< detailed description  */
				double MaxEta;  /*!< detailed description  */
				double MinR;    /*!< detailed description  */
				double MaxR;    /*!< detailed description  */

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
				double MinEt;  /*!< detailed description  */

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
