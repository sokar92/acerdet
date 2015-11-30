//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul  1 17:25:59 2015 by ROOT version 5.34/30
// from TTree ACDTree/ACDTree
// found on file: new_a-det/pythia.pp.to.Z.to.ee.conf.root
//////////////////////////////////////////////////////////

#ifndef Z_ee_h
#define Z_ee_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class Z_ee {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           ProcessID;
   Int_t           part_n;
   vector<int>     *part_pdgId;
   vector<int>     *part_mother_pdgId;
   vector<float>   *part_px;
   vector<float>   *part_py;
   vector<float>   *part_pz;
   vector<float>   *part_E;
   Int_t           ele_n;
   vector<int>     *ele_pdgId;
   vector<float>   *ele_px;
   vector<float>   *ele_py;
   vector<float>   *ele_pz;
   vector<float>   *ele_E;
   Int_t           muo_n;
   vector<int>     *muo_pdgId;
   vector<float>   *muo_px;
   vector<float>   *muo_py;
   vector<float>   *muo_pz;
   vector<float>   *muo_E;
   Int_t           pho_n;
   vector<int>     *pho_pdgId;
   vector<float>   *pho_px;
   vector<float>   *pho_py;
   vector<float>   *pho_pz;
   vector<float>   *pho_E;
   Int_t           jet_n;
   vector<int>     *jet_pdgId;
   vector<float>   *jet_px;
   vector<float>   *jet_py;
   vector<float>   *jet_pz;
   vector<float>   *jet_E;
   Float_t         pxmiss;
   Float_t         pymiss;
   Float_t         pxnue;
   Float_t         pynue;
   Float_t         pxcalo;
   Float_t         pycalo;

   // List of branches
   TBranch        *b_ProcessID;   //!
   TBranch        *b_part_n;   //!
   TBranch        *b_part_pdgId;   //!
   TBranch        *b_part_mother_pdgId;   //!
   TBranch        *b_part_px;   //!
   TBranch        *b_part_py;   //!
   TBranch        *b_part_pz;   //!
   TBranch        *b_part_E;   //!
   TBranch        *b_ele_n;   //!
   TBranch        *b_ele_pdgId;   //!
   TBranch        *b_ele_px;   //!
   TBranch        *b_ele_py;   //!
   TBranch        *b_ele_pz;   //!
   TBranch        *b_ele_E;   //!
   TBranch        *b_muo_n;   //!
   TBranch        *b_muo_pdgId;   //!
   TBranch        *b_muo_px;   //!
   TBranch        *b_muo_py;   //!
   TBranch        *b_muo_pz;   //!
   TBranch        *b_muo_E;   //!
   TBranch        *b_pho_n;   //!
   TBranch        *b_pho_pdgId;   //!
   TBranch        *b_pho_px;   //!
   TBranch        *b_pho_py;   //!
   TBranch        *b_pho_pz;   //!
   TBranch        *b_pho_E;   //!
   TBranch        *b_jet_n;   //!
   TBranch        *b_jet_pdgId;   //!
   TBranch        *b_jet_px;   //!
   TBranch        *b_jet_py;   //!
   TBranch        *b_jet_pz;   //!
   TBranch        *b_jet_E;   //!
   TBranch        *b_pxmiss;   //!
   TBranch        *b_pymiss;   //!
   TBranch        *b_pxnue;   //!
   TBranch        *b_pynue;   //!
   TBranch        *b_pxcalo;   //!
   TBranch        *b_pycalo;   //!

   Z_ee(TTree *tree=0);
   virtual ~Z_ee();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Z_ee_cxx
Z_ee::Z_ee(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../conf/pythia.pp.to.Z.to.ee.conf.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../conf/pythia.pp.to.Z.to.ee.conf.root");
      }
      f->GetObject("ACDTree",tree);

   }
   Init(tree);
}

Z_ee::~Z_ee()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Z_ee::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Z_ee::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Z_ee::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   part_pdgId = 0;
   part_mother_pdgId = 0;
   part_px = 0;
   part_py = 0;
   part_pz = 0;
   part_E = 0;
   ele_pdgId = 0;
   ele_px = 0;
   ele_py = 0;
   ele_pz = 0;
   ele_E = 0;
   muo_pdgId = 0;
   muo_px = 0;
   muo_py = 0;
   muo_pz = 0;
   muo_E = 0;
   pho_pdgId = 0;
   pho_px = 0;
   pho_py = 0;
   pho_pz = 0;
   pho_E = 0;
   jet_pdgId = 0;
   jet_px = 0;
   jet_py = 0;
   jet_pz = 0;
   jet_E = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ProcessID", &ProcessID, &b_ProcessID);
   fChain->SetBranchAddress("part_n", &part_n, &b_part_n);
   fChain->SetBranchAddress("part_pdgId", &part_pdgId, &b_part_pdgId);
   fChain->SetBranchAddress("part_mother_pdgId", &part_mother_pdgId, &b_part_mother_pdgId);
   fChain->SetBranchAddress("part_px", &part_px, &b_part_px);
   fChain->SetBranchAddress("part_py", &part_py, &b_part_py);
   fChain->SetBranchAddress("part_pz", &part_pz, &b_part_pz);
   fChain->SetBranchAddress("part_E", &part_E, &b_part_E);
   fChain->SetBranchAddress("ele_n", &ele_n, &b_ele_n);
   fChain->SetBranchAddress("ele_pdgId", &ele_pdgId, &b_ele_pdgId);
   fChain->SetBranchAddress("ele_px", &ele_px, &b_ele_px);
   fChain->SetBranchAddress("ele_py", &ele_py, &b_ele_py);
   fChain->SetBranchAddress("ele_pz", &ele_pz, &b_ele_pz);
   fChain->SetBranchAddress("ele_E", &ele_E, &b_ele_E);
   fChain->SetBranchAddress("muo_n", &muo_n, &b_muo_n);
   fChain->SetBranchAddress("muo_pdgId", &muo_pdgId, &b_muo_pdgId);
   fChain->SetBranchAddress("muo_px", &muo_px, &b_muo_px);
   fChain->SetBranchAddress("muo_py", &muo_py, &b_muo_py);
   fChain->SetBranchAddress("muo_pz", &muo_pz, &b_muo_pz);
   fChain->SetBranchAddress("muo_E", &muo_E, &b_muo_E);
   fChain->SetBranchAddress("pho_n", &pho_n, &b_pho_n);
   fChain->SetBranchAddress("pho_pdgId", &pho_pdgId, &b_pho_pdgId);
   fChain->SetBranchAddress("pho_px", &pho_px, &b_pho_px);
   fChain->SetBranchAddress("pho_py", &pho_py, &b_pho_py);
   fChain->SetBranchAddress("pho_pz", &pho_pz, &b_pho_pz);
   fChain->SetBranchAddress("pho_E", &pho_E, &b_pho_E);
   fChain->SetBranchAddress("jet_n", &jet_n, &b_jet_n);
   fChain->SetBranchAddress("jet_pdgId", &jet_pdgId, &b_jet_pdgId);
   fChain->SetBranchAddress("jet_px", &jet_px, &b_jet_px);
   fChain->SetBranchAddress("jet_py", &jet_py, &b_jet_py);
   fChain->SetBranchAddress("jet_pz", &jet_pz, &b_jet_pz);
   fChain->SetBranchAddress("jet_E", &jet_E, &b_jet_E);
   fChain->SetBranchAddress("pxmiss", &pxmiss, &b_pxmiss);
   fChain->SetBranchAddress("pymiss", &pymiss, &b_pymiss);
   fChain->SetBranchAddress("pxnue", &pxnue, &b_pxnue);
   fChain->SetBranchAddress("pynue", &pynue, &b_pynue);
   fChain->SetBranchAddress("pxcalo", &pxcalo, &b_pxcalo);
   fChain->SetBranchAddress("pycalo", &pycalo, &b_pycalo);
   Notify();
}

Bool_t Z_ee::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Z_ee::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Z_ee::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Z_ee_cxx
