#include "Configuration.h"
using namespace AcerDet::conf;

#include <fstream>
#include <sstream>

#include <vector>
#include <algorithm>
using namespace std;

/*
  Helper functions
*/
struct _split {
	string prim, sec, tail;
	_split( const string& p, const string& s, const string& t) : prim(p), sec(s), tail(t) {}
};

inline bool checkStringContains( const string& text, char c ) {
	for (string::const_iterator it = text.begin(); it != text.end(); it++)
		if ((*it) == c) return true;
	return false;
}

pair<bool,_split> split( const string& text ) {
	if (!checkStringContains(text, '.'))
		return make_pair(false, _split("","",""));
	
	char t1[128], t2[128];
	sscanf(text.c_str(), "%[^.].%[^\n]", t1, t2);
	if (!checkStringContains(string(t2), ' '))
		return make_pair(false, _split("","",""));
	
	char q1[128], q2[128];
	sscanf(t2, "%[^ ] %[^\n]", q1, q2);

	return make_pair(true, _split(string(t1), string(q1), string(q2)));
}

/*
  Following constants are binding between c++ code and *.dat configuration
*/

// -------------------------------------------------------
// -- Flag -----------------------------------------------
// -------------------------------------------------------

const std::string C_Flag					= "Flag";
const std::string C_Flag_HistogramID		= "HistogramID";
const std::string C_Flag_Smearing			= "Smearing";
const std::string C_Flag_BField				= "B-Field";
const std::string C_Flag_Susy_Lsp_Particle	= "SUSY_LSP_Particle";
const std::string C_Flag_BC_Jets_Labeling	= "BC-JetsLabeling";
const std::string C_Flag_Tau_Jets_Labeling	= "Tau-JetsLabeling";
const std::string C_Flag_JetCalibration		= "JetCalibration";

Configuration::_Flag::_Flag() :
	HistogramId( 10000 ),
	Smearing( true ),
	BField( true ),
	SusyParticle( 66 ),
	BCJetsLabeling( true ),
	TauJetsLabeling( true ),
	JetCalibration( true ) 
{}

void Configuration::_Flag::read( const string& line ) {
	pair<bool,_split> P = split(line);
	
	if (!P.first || P.second.prim != C_Flag)
		return;

	if (P.second.sec == C_Flag_HistogramID) 			sscanf(P.second.tail.c_str(), "%d", &HistogramId);
	else if (P.second.sec == C_Flag_Susy_Lsp_Particle)	sscanf(P.second.tail.c_str(), "%d", &SusyParticle);
	else if (P.second.sec == C_Flag_Smearing) 			Smearing = (P.second.tail == "1");
	else if (P.second.sec == C_Flag_BField) 			BField = (P.second.tail == "1");
	else if (P.second.sec == C_Flag_BC_Jets_Labeling)	BCJetsLabeling = (P.second.tail == "1");
	else if (P.second.sec == C_Flag_Tau_Jets_Labeling)	TauJetsLabeling = (P.second.tail == "1");
	else if (P.second.sec == C_Flag_JetCalibration)		JetCalibration = (P.second.tail == "1");
}

string Configuration::_Flag::write( ) const {
	stringstream ss;
	ss << C_Flag << "." << C_Flag_HistogramID	<< " " << HistogramId		<< endl;
	ss << C_Flag << "." << C_Flag_Smearing		<< " " << Smearing		<< endl;
	ss << C_Flag << "." << C_Flag_BField		<< " " << BField		<< endl;
	ss << C_Flag << "." << C_Flag_Susy_Lsp_Particle	<< " " << SusyParticle 		<< endl;
	ss << C_Flag << "." << C_Flag_BC_Jets_Labeling	<< " " << BCJetsLabeling 	<< endl;
	ss << C_Flag << "." << C_Flag_Tau_Jets_Labeling	<< " " << TauJetsLabeling 	<< endl;
	ss << C_Flag << "." << C_Flag_JetCalibration	<< " " << JetCalibration	<< endl;
	return ss.str();
}

// -------------------------------------------------------
// -- Cell -----------------------------------------------
// -------------------------------------------------------

const std::string C_Cell					= "Cell";
const std::string C_Cell_RapidityCoverage	= "RapidityCoverage";
const std::string C_Cell_Min_pT				= "MinpT";
const std::string C_Cell_Min_Et				= "MinEt";
const std::string C_Cell_EtaTransition		= "EtaTransition";
const std::string C_Cell_GranularityEta		= "GranularityEta";
const std::string C_Cell_GranularityPhi		= "GranularityPhi";

Configuration::_Cell::_Cell() :
	RapidityCoverage( 5.0 ),
	MinpT( 0.5 ),
	MinEt( 0.0 ),
	EtaTransition( 3.2 ),
	GranularityEta( 0.1 ),
	GranularityPhi( 0.1 )
{}

void Configuration::_Cell::read( const string& line ) {
	pair<bool,_split> P = split(line);
	
	if (!P.first || P.second.prim != C_Cell)
		return;

	if (P.second.sec == C_Cell_RapidityCoverage) 	sscanf(P.second.tail.c_str(), "%lf", &RapidityCoverage);
	else if (P.second.sec == C_Cell_Min_pT)			sscanf(P.second.tail.c_str(), "%lf", &MinpT);
	else if (P.second.sec == C_Cell_Min_Et)			sscanf(P.second.tail.c_str(), "%lf", &MinEt);
	else if (P.second.sec == C_Cell_EtaTransition)	sscanf(P.second.tail.c_str(), "%lf", &EtaTransition);
	else if (P.second.sec == C_Cell_GranularityEta)	sscanf(P.second.tail.c_str(), "%lf", &GranularityEta);
	else if (P.second.sec == C_Cell_GranularityPhi)	sscanf(P.second.tail.c_str(), "%lf", &GranularityPhi);
}

string Configuration::_Cell::write( ) const {
	stringstream ss;
	ss << C_Cell << "." << C_Cell_RapidityCoverage	<< " " << RapidityCoverage	<< endl;
	ss << C_Cell << "." << C_Cell_Min_pT		<< " " << MinpT			<< endl;
	ss << C_Cell << "." << C_Cell_Min_Et		<< " " << MinEt			<< endl;
	ss << C_Cell << "." << C_Cell_EtaTransition	<< " " << EtaTransition		<< endl;
	ss << C_Cell << "." << C_Cell_GranularityEta	<< " " << GranularityEta	<< endl;
	ss << C_Cell << "." << C_Cell_GranularityPhi	<< " " << GranularityPhi	<< endl;
	return ss.str();
}

// -------------------------------------------------------
// -- Cluster --------------------------------------------
// -------------------------------------------------------

const std::string C_Cluster					= "Cluster";
const std::string C_Cluster_Min_Et			= "MinEt";
const std::string C_Cluster_ConeR			= "ConeR";
const std::string C_Cluster_RapidityCoverage= "RapidityCoverage";
const std::string C_Cluster_Min_Et_init		= "MinEtInit";

Configuration::_Cluster::_Cluster() :
	RapidityCoverage( 5.0 ),
	ConeR( 0.4 ),
	MinEt( 5.0 ),
	MinEtInit( 1.5 )
{}

void Configuration::_Cluster::read( const string& line ) {
	pair<bool,_split> P = split(line);
	
	if (!P.first || P.second.prim != C_Cluster)
		return;

	if (P.second.sec == C_Cluster_Min_Et) 				sscanf(P.second.tail.c_str(), "%lf", &MinEt);
	else if (P.second.sec == C_Cluster_ConeR)			sscanf(P.second.tail.c_str(), "%lf", &ConeR);
	else if (P.second.sec == C_Cluster_RapidityCoverage)sscanf(P.second.tail.c_str(), "%lf", &RapidityCoverage);
	else if (P.second.sec == C_Cluster_Min_Et_init)		sscanf(P.second.tail.c_str(), "%lf", &MinEtInit);
}

string Configuration::_Cluster::write( ) const {
	stringstream ss;
	ss << C_Cluster << "." << C_Cluster_Min_Et		<< " " << MinEt			<< endl;
	ss << C_Cluster << "." << C_Cluster_ConeR		<< " " << ConeR			<< endl;
	ss << C_Cluster << "." << C_Cluster_RapidityCoverage	<< " " << RapidityCoverage	<< endl;
	ss << C_Cluster << "." << C_Cluster_Min_Et_init		<< " " << MinEtInit		<< endl;
	return ss.str();
}

// -------------------------------------------------------
// -- Muon -----------------------------------------------
// -------------------------------------------------------

const std::string C_Muon				= "Muon";
const std::string C_Muon_MinMomenta		= "MinMomenta";
const std::string C_Muon_MaxEta			= "MaxEta";
const std::string C_Muon_MinIsolRlj		= "MinIsolRlj";
const std::string C_Muon_ConeR			= "ConeR";
const std::string C_Muon_MaxEnergy		= "MaxEnergy";

Configuration::_Muon::_Muon() :
	MinMomenta( 6.0 ),
	MaxEta( 2.5 ),
	MinIsolRlj( 0.4 ),
	ConeR( 0.2 ),
	MaxEnergy( 10.0 )
{}

void Configuration::_Muon::read( const string& line ) {
	pair<bool,_split> P = split(line);
	
	if (!P.first || P.second.prim != C_Muon)
		return;

	if (P.second.sec == C_Muon_MinMomenta) 		sscanf(P.second.tail.c_str(), "%lf", &MinMomenta);
	else if (P.second.sec == C_Muon_MaxEta)		sscanf(P.second.tail.c_str(), "%lf", &MaxEta);
	else if (P.second.sec == C_Muon_MinIsolRlj) sscanf(P.second.tail.c_str(), "%lf", &MinIsolRlj);
	else if (P.second.sec == C_Muon_ConeR)		sscanf(P.second.tail.c_str(), "%lf", &ConeR);
	else if (P.second.sec == C_Muon_MaxEnergy)	sscanf(P.second.tail.c_str(), "%lf", &MaxEnergy);
}

string Configuration::_Muon::write( ) const {
	stringstream ss;
	ss << C_Muon << "." << C_Muon_MinMomenta	<< " " << MinMomenta	<< endl;
	ss << C_Muon << "." << C_Muon_MaxEta		<< " " << MaxEta	<< endl;
	ss << C_Muon << "." << C_Muon_MinIsolRlj	<< " " << MinIsolRlj	<< endl;
	ss << C_Muon << "." << C_Muon_ConeR		<< " " << ConeR		<< endl;
	ss << C_Muon << "." << C_Muon_MaxEnergy		<< " " << MaxEnergy	<< endl;
	return ss.str();
}

// -------------------------------------------------------
// -- Photon ---------------------------------------------
// -------------------------------------------------------

const std::string C_Photon				= "Photon";
const std::string C_Photon_MinMomenta	= "MinMomenta";
const std::string C_Photon_MaxEta		= "MaxEta";
const std::string C_Photon_MinJetsRlj	= "MinJetsRlj";
const std::string C_Photon_MinIsolRlj	= "MinIsolRlj";
const std::string C_Photon_ConeR		= "ConeR";
const std::string C_Photon_MaxEnergy	= "MaxEnergy";

Configuration::_Photon::_Photon() :
	MinMomenta( 5.0 ),
	MaxEta( 2.5 ),
	MinJetsRlj( 0.15 ),
	MinIsolRlj( 0.4 ),
	ConeR( 0.2 ),
	MaxEnergy( 10.0 )
{}

void Configuration::_Photon::read( const string& line ) {
	pair<bool,_split> P = split(line);
	
	if (!P.first || P.second.prim != C_Photon)
		return;

	if (P.second.sec == C_Photon_MinMomenta) 		sscanf(P.second.tail.c_str(), "%lf", &MinMomenta);
	else if (P.second.sec == C_Photon_MaxEta)		sscanf(P.second.tail.c_str(), "%lf", &MaxEta);
	else if (P.second.sec == C_Photon_MinJetsRlj) 	sscanf(P.second.tail.c_str(), "%lf", &MinJetsRlj);
	else if (P.second.sec == C_Photon_MinIsolRlj) 	sscanf(P.second.tail.c_str(), "%lf", &MinIsolRlj);
	else if (P.second.sec == C_Photon_ConeR)		sscanf(P.second.tail.c_str(), "%lf", &ConeR);
	else if (P.second.sec == C_Photon_MaxEnergy)	sscanf(P.second.tail.c_str(), "%lf", &MaxEnergy);
}

string Configuration::_Photon::write( ) const {
	stringstream ss;
	ss << C_Photon << "." << C_Photon_MinMomenta	<< " " << MinMomenta	<< endl;
	ss << C_Photon << "." << C_Photon_MaxEta	<< " " << MaxEta	<< endl;
	ss << C_Photon << "." << C_Photon_MinJetsRlj	<< " " << MinJetsRlj	<< endl;
	ss << C_Photon << "." << C_Photon_MinIsolRlj	<< " " << MinIsolRlj	<< endl;
	ss << C_Photon << "." << C_Photon_ConeR		<< " " << ConeR		<< endl;
	ss << C_Photon << "." << C_Photon_MaxEnergy	<< " " << MaxEnergy	<< endl;
	return ss.str();
}

// -------------------------------------------------------
// -- Electron -------------------------------------------
// -------------------------------------------------------

const std::string C_Electron			= "Electron";
const std::string C_Electron_MinMomenta	= "MinMomenta";
const std::string C_Electron_MaxEta		= "MaxEta";
const std::string C_Electron_MinJetsRlj	= "MinJetsRlj";
const std::string C_Electron_MinIsolRlj	= "MinIsolRlj";
const std::string C_Electron_ConeR		= "ConeR";
const std::string C_Electron_MaxEnergy	= "MaxEnergy";

Configuration::_Electron::_Electron() :
	MinMomenta( 5.0 ),
	MaxEta( 2.5 ),
	MinJetsRlj( 0.15 ),
	MinIsolRlj( 0.4 ),
	ConeR( 0.2 ),
	MaxEnergy( 10.0 )
{}

void Configuration::_Electron::read( const string& line ) {
	pair<bool,_split> P = split(line);
	
	if (!P.first || P.second.prim != C_Electron)
		return;

	if (P.second.sec == C_Electron_MinMomenta) 		sscanf(P.second.tail.c_str(), "%lf", &MinMomenta);
	else if (P.second.sec == C_Electron_MaxEta)		sscanf(P.second.tail.c_str(), "%lf", &MaxEta);
	else if (P.second.sec == C_Electron_MinJetsRlj) sscanf(P.second.tail.c_str(), "%lf", &MinJetsRlj);
	else if (P.second.sec == C_Electron_MinIsolRlj) sscanf(P.second.tail.c_str(), "%lf", &MinIsolRlj);
	else if (P.second.sec == C_Electron_ConeR)		sscanf(P.second.tail.c_str(), "%lf", &ConeR);
	else if (P.second.sec == C_Electron_MaxEnergy)	sscanf(P.second.tail.c_str(), "%lf", &MaxEnergy);
}

string Configuration::_Electron::write( ) const {
	stringstream ss;
	ss << C_Electron << "." << C_Electron_MinMomenta	<< " " << MinMomenta	<< endl;
	ss << C_Electron << "." << C_Electron_MaxEta		<< " " << MaxEta	<< endl;
	ss << C_Electron << "." << C_Electron_MinJetsRlj	<< " " << MinJetsRlj	<< endl;
	ss << C_Electron << "." << C_Electron_MinIsolRlj	<< " " << MinIsolRlj	<< endl;
	ss << C_Electron << "." << C_Electron_ConeR		<< " " << ConeR		<< endl;
	ss << C_Electron << "." << C_Electron_MaxEnergy		<< " " << MaxEnergy	<< endl;
	return ss.str();
}

// -------------------------------------------------------
// -- Jet ------------------------------------------------
// ------------------------------------------------------- 

const std::string C_Jets					= "Jets";
const std::string C_Jets_MinEnergy			= "MinEnergy";
const std::string C_Jets_RapidityCoverage	= "RapidityCoverage";

Configuration::_Jet::_Jet() :
	RapidityCoverage( 5.0 ),
	MinEnergy( 10.0 )
{}

void Configuration::_Jet::read( const string& line ) {
	pair<bool,_split> P = split(line);
	
	if (!P.first || P.second.prim != C_Jets)
		return;

	if (P.second.sec == C_Jets_MinEnergy)	 			sscanf(P.second.tail.c_str(), "%lf", &MinEnergy);
	else if (P.second.sec == C_Jets_RapidityCoverage)	sscanf(P.second.tail.c_str(), "%lf", &RapidityCoverage);
}

string Configuration::_Jet::write( ) const {
	stringstream ss;
	ss << C_Jets << "." << C_Jets_MinEnergy		<< " " << MinEnergy		<< endl;
	ss << C_Jets << "." << C_Jets_RapidityCoverage	<< " " << RapidityCoverage	<< endl;
	return ss.str();
}

// -------------------------------------------------------
// -- BJet -----------------------------------------------
// -------------------------------------------------------

const std::string C_BJets			= "BJets";
const std::string C_BJets_MinMomenta= "MinMomenta";
const std::string C_BJets_MaxEta	= "MaxEta";
const std::string C_BJets_MaxRbj	= "MaxRbj";

Configuration::_BJet::_BJet() :
	MinMomenta( 5.0 ),
	MaxEta( 2.5 ),
	MaxRbj( 0.2 )
{}

void Configuration::_BJet::read( const string& line ) {
	pair<bool,_split> P = split(line);
	
	if (!P.first || P.second.prim != C_BJets)
		return;
	
	if (P.second.sec == C_BJets_MinMomenta) 	sscanf(P.second.tail.c_str(), "%lf", &MinMomenta);
	else if (P.second.sec == C_BJets_MaxEta)	sscanf(P.second.tail.c_str(), "%lf", &MaxEta);
	else if (P.second.sec == C_BJets_MaxRbj)	sscanf(P.second.tail.c_str(), "%lf", &MaxRbj);
}

string Configuration::_BJet::write( ) const {
	stringstream ss;
	ss << C_BJets << "." << C_BJets_MinMomenta	<< " " << MinMomenta	<< endl;
	ss << C_BJets << "." << C_BJets_MaxEta		<< " " << MaxEta	<< endl;
	ss << C_BJets << "." << C_BJets_MaxRbj		<< " " << MaxRbj	<< endl;
	return ss.str();
}

// -------------------------------------------------------
// -- CJets ----------------------------------------------
// -------------------------------------------------------

const std::string C_CJets			= "CJets";
const std::string C_CJets_MinMomenta= "MinMomenta";
const std::string C_CJets_MaxEta	= "MaxEta";
const std::string C_CJets_MaxRcj	= "MaxRcj";

Configuration::_CJet::_CJet() :
	MinMomenta( 5.0 ),
	MaxEta( 2.5 ),
	MaxRcj( 0.2 )
{}

void Configuration::_CJet::read( const string& line ) {
	pair<bool,_split> P = split(line);
	
	if (!P.first || P.second.prim != C_CJets)
		return;

	if (P.second.sec == C_CJets_MinMomenta) 	sscanf(P.second.tail.c_str(), "%lf", &MinMomenta);
	else if (P.second.sec == C_CJets_MaxEta)	sscanf(P.second.tail.c_str(), "%lf", &MaxEta);
	else if (P.second.sec == C_CJets_MaxRcj)	sscanf(P.second.tail.c_str(), "%lf", &MaxRcj);
}

string Configuration::_CJet::write( ) const {
	stringstream ss;
	ss << C_CJets << "." << C_CJets_MinMomenta	<< " " << MinMomenta	<< endl;
	ss << C_CJets << "." << C_CJets_MaxEta		<< " " << MaxEta	<< endl;
	ss << C_CJets << "." << C_CJets_MaxRcj		<< " " << MaxRcj	<< endl;
	return ss.str();
}

// -------------------------------------------------------
// -- Tau ------------------------------------------------
// -------------------------------------------------------

const std::string C_Tau				= "Tau";
const std::string C_Tau_MinpT		= "MinpT";
const std::string C_Tau_MaxEta		= "MaxEta";
const std::string C_Tau_MinR		= "MinR";
const std::string C_Tau_MaxR		= "MaxR";

Configuration::_Tau::_Tau() :
	MinpT( 10.0 ),
	MaxEta( 2.5 ),
	MinR( 0.3 ),
	MaxR( 0.9 )
{}

void Configuration::_Tau::read( const string& line ) {
	pair<bool,_split> P = split(line);
	
	if (!P.first || P.second.prim != C_Tau)
		return;

	if (P.second.sec == C_Tau_MinpT)	 	sscanf(P.second.tail.c_str(), "%lf", &MinpT);
	else if (P.second.sec == C_Tau_MaxEta)	sscanf(P.second.tail.c_str(), "%lf", &MaxEta);
	else if (P.second.sec == C_Tau_MinR)	sscanf(P.second.tail.c_str(), "%lf", &MinR);
	else if (P.second.sec == C_Tau_MaxR)	sscanf(P.second.tail.c_str(), "%lf", &MaxR);
}

string Configuration::_Tau::write( ) const {
	stringstream ss;
	ss << C_Tau << "." << C_Tau_MinpT	<< " " << MinpT		<< endl;
	ss << C_Tau << "." << C_Tau_MaxEta	<< " " << MaxEta	<< endl;
	ss << C_Tau << "." << C_Tau_MinR	<< " " << MinR		<< endl;
	ss << C_Tau << "." << C_Tau_MaxR	<< " " << MaxR		<< endl;
	return ss.str();
}

// -------------------------------------------------------
// -- Misc -----------------------------------------------
// -------------------------------------------------------

const std::string C_Misc			= "Misc";
const std::string C_Misc_MinEt		= "MinEt";

Configuration::_Misc::_Misc() :
	MinEt( 0.0 )
{}

void Configuration::_Misc::read( const string& line ) {
	pair<bool,_split> P = split(line);
	
	if (!P.first || P.second.prim != C_Misc)
		return;

	if (P.second.sec == C_Misc_MinEt) 
		sscanf(P.second.tail.c_str(), "%lf", &MinEt);
}

string Configuration::_Misc::write( ) const {
	stringstream ss;
	ss << C_Misc << "." << C_Misc_MinEt	<< " " << MinEt	<< endl;
	return ss.str();
}


// -------------------------------------------------------
// -- AcerDetConfiguration members -----------------------
// -------------------------------------------------------


/*
  Create default configuration
*/
Configuration Configuration::getDefault() {
	Configuration config;
	return config;
}

/*
  Read configuration from *.dat file
*/
Configuration Configuration::fromFile( const std::string& fileName ) {
	Configuration config;

	fstream stream;
	stream.open(fileName.c_str(), ios_base::in);
	if (stream.is_open()) {
		
		string line;
		while (stream) {
			std::getline(stream, line);

			config.Flag.read( line );
			config.Cell.read( line );
			config.Cluster.read( line );
			config.Muon.read( line );
			config.Photon.read( line );
			config.Electron.read( line );
			config.Jet.read( line );
			config.BJet.read( line );
			config.CJet.read( line );
			config.Tau.read( line );
			config.Misc.read( line );
		}

		stream.close();
	} else {
		printf ("Acer DET 2.0 Configuration read -> open file failed!\n");
		throw -1;
	}

	return config;
}

/*
  Save configuration into file
*/
void Configuration::save( const Configuration& config, const std::string& fileName ) {
	std::fstream stream;
	stream.open(fileName.c_str(), ios_base::out);
	if (stream.is_open()) {
		stream 	<< config.Flag.write() 		<< endl
			<< config.Cell.write() 		<< endl
			<< config.Cluster.write()	<< endl
			<< config.Muon.write()		<< endl
			<< config.Photon.write()	<< endl
			<< config.Electron.write()	<< endl
			<< config.Jet.write()		<< endl
			<< config.BJet.write()		<< endl
			<< config.CJet.write()		<< endl
			<< config.Tau.write()		<< endl
			<< config.Misc.write();
		stream.close();
	} else {
		printf ("Acer DET 2.0 Configuration save -> open file failed\n");
		throw -1;
	}
}
