#include "Root_NTupleManager.h"
#include <cstring>
using namespace AcerDet::external;

#include "../src/core/Functions.h"
using namespace AcerDet::core;

Root_NTupleManager::Root_NTupleManager()
{}

Root_NTupleManager::~Root_NTupleManager() 
{}

void Root_NTupleManager::init() {
	
	/* in order to enable using std::vectors */
	gROOT->ProcessLine("#include<vector>");
	
	ntuple = new TTree("ACDTree", "ACDTree");
	
	ntuple->Branch("ProcessID",          &n_ProcessID);
	
	ntuple->Branch("part_n",             &n_part_n);
	ntuple->Branch("part_pdgId",         &n_part_pdgId);
	ntuple->Branch("part_px",            &n_part_px);
	ntuple->Branch("part_py",            &n_part_py);
	ntuple->Branch("part_pz",            &n_part_pz);
	ntuple->Branch("part_E",             &n_part_E);

	ntuple->Branch("ele_n",              &n_ele_n);
	ntuple->Branch("ele_pdgId",          &n_ele_pdgId);
	ntuple->Branch("ele_px",             &n_ele_px);
	ntuple->Branch("ele_py",             &n_ele_py);
	ntuple->Branch("ele_pz",             &n_ele_pz);
	ntuple->Branch("ele_E",              &n_ele_E);

	ntuple->Branch("muo_n",              &n_muo_n);
	ntuple->Branch("muo_pdgId",          &n_muo_pdgId);
	ntuple->Branch("muo_px",             &n_muo_px);
	ntuple->Branch("muo_py",             &n_muo_py);
	ntuple->Branch("muo_pz",             &n_muo_pz);
	ntuple->Branch("muo_E",              &n_muo_E);

	ntuple->Branch("pho_n",              &n_pho_n);
	ntuple->Branch("pho_pdgId",          &n_pho_pdgId);
	ntuple->Branch("pho_px",             &n_pho_px);
	ntuple->Branch("pho_py",             &n_pho_py);
	ntuple->Branch("pho_pz",             &n_pho_pz);
	ntuple->Branch("pho_E",              &n_pho_E);

	ntuple->Branch("jet_n",              &n_jet_n);
	ntuple->Branch("jet_pdgId",          &n_jet_pdgId);
	ntuple->Branch("jet_px",             &n_jet_px);
	ntuple->Branch("jet_py",             &n_jet_py);
	ntuple->Branch("jet_pz",             &n_jet_pz);
	ntuple->Branch("jet_E",              &n_jet_E);

	ntuple->Branch("pxmiss",             &n_pxmiss);
	ntuple->Branch("pymiss",             &n_pxmiss);
	ntuple->Branch("pxnue",              &n_pxnue);
	ntuple->Branch("pynue",              &n_pynue);
	ntuple->Branch("pxcalo",             &n_pxcalo);
	ntuple->Branch("pycalo",             &n_pycalo);
}

void Root_NTupleManager::fill(
	const InputRecord& irecord,
	const OutputRecord& orecord,
	Real64_t weigth
) {
	n_part_px.clear();
	n_part_py.clear();
	n_part_pz.clear();
	n_part_E.clear();
	n_part_pdgId.clear();
	
	int count = 0;
	const vector<Particle>& parts = irecord.particles();
	for (int j=0; j<parts.size(); ++j) {
		n_part_pdgId.push_back(parts[j].pdg_id);
		
		// pick only particles from hard process
		if (!isHardProcess(parts, j))
			continue;
			
		Real64_t px = parts[j].pT() * cos(parts[j].getPhi());
		Real64_t py = parts[j].pT() * sin(parts[j].getPhi());
		Real64_t pz = parts[j].pT() * sinh(parts[j].getEta());
		Real64_t E  = parts[j].pT() * cosh(parts[j].getEta());
		
		n_part_px.push_back(px);
		n_part_py.push_back(py);
		n_part_pz.push_back(pz);
		n_part_E.push_back(E);
		count++;
	}
	n_part_n = count;   

	n_pho_px.clear();
	n_pho_py.clear();
	n_pho_pz.clear();
	n_pho_E.clear();
	n_pho_pdgId.clear();
	
	for (int j=0; j<orecord.Photons.size(); ++j) {
		n_pho_pdgId.push_back(orecord.Photons[j].pdg_id);
		
		Real64_t px = orecord.Photons[j].pT * cos(orecord.Photons[j].phi);
		Real64_t py = orecord.Photons[j].pT * sin(orecord.Photons[j].phi);
		Real64_t pz = orecord.Photons[j].pT * sinh(orecord.Photons[j].eta);
		Real64_t E  = orecord.Photons[j].pT * cosh(orecord.Photons[j].eta);
		
		n_pho_px.push_back(px);
		n_pho_py.push_back(py);
		n_pho_pz.push_back(pz);
		n_pho_E.push_back(E);
	}
	n_pho_n = orecord.Photons.size(); 

	n_muo_px.clear();
	n_muo_py.clear();
	n_muo_pz.clear();
	n_muo_E.clear();
	n_muo_pdgId.clear();
	
	for (int j=0; j<orecord.Muons.size(); ++j) {
		n_muo_pdgId.push_back(orecord.Muons[j].pdg_id);
		
		Real64_t px = orecord.Muons[j].pT * cos(orecord.Muons[j].phi);
		Real64_t py = orecord.Muons[j].pT * sin(orecord.Muons[j].phi);
		Real64_t pz = orecord.Muons[j].pT * sinh(orecord.Muons[j].eta);
		Real64_t E  = orecord.Muons[j].pT * cosh(orecord.Muons[j].eta);
		
		n_muo_px.push_back(px);
		n_muo_py.push_back(py);
		n_muo_pz.push_back(pz);
		n_muo_E.push_back(E);
	}
	n_muo_n = orecord.Muons.size(); 

	n_ele_px.clear();
	n_ele_py.clear();
	n_ele_pz.clear();
	n_ele_E.clear();
	n_ele_pdgId.clear();
	
	for (int j=0; j<orecord.Electrons.size(); ++j) {
		n_ele_pdgId.push_back(orecord.Electrons[j].pdg_id);
		
		Real64_t px = orecord.Electrons[j].pT * cos(orecord.Electrons[j].phi);
		Real64_t py = orecord.Electrons[j].pT * sin(orecord.Electrons[j].phi);
		Real64_t pz = orecord.Electrons[j].pT * sinh(orecord.Electrons[j].eta);
		Real64_t E  = orecord.Electrons[j].pT * cosh(orecord.Electrons[j].eta);
		
		n_ele_px.push_back(px);
		n_ele_py.push_back(py);
		n_ele_pz.push_back(pz);
		n_ele_E.push_back(E);
	}
	n_ele_n = orecord.Electrons.size(); 
	
	n_jet_px.clear();
	n_jet_py.clear();
	n_jet_pz.clear();
	n_jet_E.clear();
	n_jet_pdgId.clear();
	
	for (int j=0; j<orecord.Jets.size(); ++j) {
		n_jet_pdgId.push_back(orecord.Jets[j].type);
		
		Real64_t px = orecord.Jets[j].pT * cos(orecord.Jets[j].phi);
		Real64_t py = orecord.Jets[j].pT * sin(orecord.Jets[j].phi);
		Real64_t pz = orecord.Jets[j].pT * sinh(orecord.Jets[j].eta);
		Real64_t E  = orecord.Jets[j].pT * cosh(orecord.Jets[j].eta);
		
		n_jet_px.push_back(px);
		n_jet_py.push_back(py);
		n_jet_pz.push_back(pz);
		n_jet_E.push_back(E);
	}
	n_jet_n = orecord.Jets.size(); 

	ntuple->Fill();
}

void Root_NTupleManager::write() {
	ntuple->Write();
}
