#define enu_cxx
#include "W_enu.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void enu::Loop()
{

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   TH1F *px_miss = new TH1F("shist12602" , "Mis: Px_miss - Px_nu", 100, -100, 100 );
   TH1F *py_miss = new TH1F("shist12607" , "Mis: Py_miss - Py_nu", 100, -100, 100 );
   TH1F *px_calo = new TH1F("shist12612" , "Mis: Px_calo - Px_nu", 100, -100, 100 );
   TH1F *py_calo = new TH1F("shist12617" , "Mis: Py_calo - Py_nu", 100, -100, 100 );
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      GetEntry(jentry);
      px_miss -> Fill ( pxmiss - pxnue );
      py_miss -> Fill ( pymiss - pynue );
      px_calo -> Fill ( pxcalo - pxnue );
      py_calo -> Fill ( pycalo - pynue );

   }
}
