#ifndef _AcerDet_External_Root_NTupleManager_
#define _AcerDet_External_Root_NTupleManager_

#include <TROOT.h>
#include <TTree.h>

#include "../src/io/OutputRecord.h"
#include "../src/io/InputRecord.h"
using namespace AcerDet::io;

namespace AcerDet {
	namespace external {
		
		class Root_NTupleManager {
		private:
			//ACDtree
			TTree* ntuple;
			
			Int32_t n_part_n;
			vector<Real32_t> n_part_px;
			vector<Real32_t> n_part_py;
			vector<Real32_t> n_part_pz;
			vector<Real32_t> n_part_E;
			vector<Int32_t> n_part_pdgId;

			Int32_t  n_muo_n;
			vector<Real32_t> n_muo_px;
			vector<Real32_t> n_muo_py;
			vector<Real32_t> n_muo_pz;
			vector<Real32_t> n_muo_E;
			vector<Int32_t> n_muo_pdgId;

			Int32_t  n_pho_n;
			vector<Real32_t> n_pho_px;
			vector<Real32_t> n_pho_py;
			vector<Real32_t> n_pho_pz;
			vector<Real32_t> n_pho_E;
			vector<Int32_t> n_pho_pdgId;
 
			Int32_t  n_ele_n;
			vector<Real32_t> n_ele_px;
			vector<Real32_t> n_ele_py;
			vector<Real32_t> n_ele_pz;
			vector<Real32_t> n_ele_E;
			vector<Int32_t> n_ele_pdgId;

			Int32_t  n_jet_n;
			vector<Real32_t> n_jet_px;
			vector<Real32_t> n_jet_py;
			vector<Real32_t> n_jet_pz;
			vector<Real32_t> n_jet_E;
			vector<Int32_t> n_jet_pdgId;

			Real32_t n_pxmiss;
			Real32_t n_pymiss;
			Real32_t n_pxnue;
			Real32_t n_pynue;
			Real32_t n_pxcalo;
			Real32_t n_pycalo;
			Int32_t n_ProcessID;
  
		public:
			Root_NTupleManager();
			
			~Root_NTupleManager();
		
			void init();
			
			void fill(
				const InputRecord& iRec,
				const OutputRecord& oRec,
				Real64_t weight );
			
			void write();
		};

	}
}

#endif
