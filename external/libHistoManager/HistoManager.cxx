#define HistManager_cxx

#include "HistoManager.h"
#include "TF1.h"
#include "TH1.h"   
#include "TH2.h"
#include "TProfile.h"

HistoManager* HistoManager::gHstMan = 0;

HistoManager* HistoManager::getInstance()
{
  if (!gHstMan) gHstMan=new HistoManager();
  return gHstMan;
}

ClassImp(HistoManager)


///////////////////////////////////////////////////////////
HistDefTH1F::HistDefTH1F(const char *name, const char *title, int bins, float xlow, float xhigh) {// constructor
  
  m_name = name;
  m_title = title;
  m_nbins = bins;
  m_xlow = xlow;
  m_xhigh = xhigh;
}

///////////////////////////////////////////////////////////
HistDefTH2F::HistDefTH2F(const char *name, const char *title, int bins1, float xlow1, float xhigh1,int bins2, float xlow2, float xhigh2) {// constructor
  
  m_name = name;
  m_title = title;
  m_n1bins = bins1;
  m_n2bins = bins2;
  m_x1low = xlow1;
  m_x1high = xhigh1;
  m_x2low = xlow2;
  m_x2high = xhigh2;
}

///////////////////////////////////////////////////////////
HistDefTH3F::HistDefTH3F(const char *name, const char *title, int bins1, float xlow1, float xhigh1,int bins2, float xlow2, float xhigh2,
                                                              int bins3, float xlow3, float xhigh3 ) {// constructor
  
  m_name = name;
  m_title = title;
  m_n1bins = bins1;
  m_n2bins = bins2;
  m_n3bins = bins3;
  m_x1low = xlow1;
  m_x1high = xhigh1;
  m_x2low = xlow2;
  m_x2high = xhigh2;
  m_x3low = xlow3;
  m_x3high = xhigh3;
}

///////////////////////////////////////////////////////////
HistDefTProf::HistDefTProf(const char *name, const char *title, int bins, float xlow1, float xhigh1,float xlow2, float xhigh2){// constructor
  
  m_name = name;
  m_title = title;
  m_nbins = bins;
  m_x1low = xlow1;
  m_x1high = xhigh1;
  m_x2low = xlow2;
  m_x2high = xhigh2;
}


//////////////////////////////////////////////////////////////////////////////////////
/// Constructor

HistoManager::HistoManager() {

  th1fSize  = 0;
  th2fSize  = 0;
  th3fSize  = 0;
  tprofSize = 0;

}

/////////////////////////////////////////////////////////////////////////////////////
/// Destructor - check up memory allocation
/// delete any memory allocation on the heap

HistoManager::~HistoManager() {

  if(m_TH1FArray)  delete m_TH1FArray;
  if(m_TH2FArray)  delete m_TH2FArray;
  if(m_TH3FArray)  delete m_TH3FArray;
  if(m_TProfArray) delete m_TProfArray;

}


/////////////////////////////////////////////////////////////////////////////////////
//
// Define  histograms for XsectionAnalysis analysis
//
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//
//Add new TH1F histogram 
//
/////////////////////////////////////////////////////////////////////////////////////
void HistoManager::addTH1F(const char *name, const char *title, int bins, float xlow, float xhigh) {
  TClonesArray &a = *m_TH1FArray;
  a[th1fSize] = new (a[th1fSize]) TH1F(name, title, bins, xlow, xhigh);
//  std::cout << "histoTH1F    " << name << std::endl;
  th1fSize++;
  
}
/////////////////////////////////////////////////////////////////////////////////////
//
//Add new TH1F histogram
//
/////////////////////////////////////////////////////////////////////////////////////
void HistoManager::addTH1F(const char *name, const char *title, int binsx, const double *xbins ) {
  TClonesArray &a = *m_TH1FArray;
  a[th1fSize] = new (a[th1fSize]) TH1F(name, title, binsx, xbins);
//  std::cout << "histoTH1F    " << name << std::endl;
  th1fSize++;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//Add new TH3F histogram
//
/////////////////////////////////////////////////////////////////////////////////////
void HistoManager::addTH3F(const char *name, const char *title, int bins1, float xlow1, float xhigh1,int bins2, float xlow2, float xhigh2,
                           int bins3, float xlow3, float xhigh3) {
  TClonesArray &a = *m_TH3FArray;
  a[th3fSize] = new (a[th3fSize]) TH3F(name, title, bins1, xlow1, xhigh1, bins2, xlow2, xhigh2, bins3, xlow3, xhigh3);
//  std::cout << "histoTH3F    " << name << std::endl;
  th3fSize++;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//Add new TH3F histogram
//
/////////////////////////////////////////////////////////////////////////////////////
void HistoManager::addTH3F(const char *name, const char *title, int binsx, const double *xbins , int binsy, const double *ybins, 
                           int binsz, const double *zbins) {
  TClonesArray &a = *m_TH3FArray;
  a[th3fSize] = new (a[th3fSize]) TH3F(name, title, binsx, xbins, binsy, ybins, binsz, zbins);
//  std::cout << "histoTH3F    " << name << std::endl;
  th3fSize++;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//Add new TH2F histogram
//
/////////////////////////////////////////////////////////////////////////////////////
void HistoManager::addTH2F(const char *name, const char *title, int bins1, float xlow1, float xhigh1,int bins2, float xlow2, float xhigh2) {
  TClonesArray &a = *m_TH2FArray;
  a[th2fSize] = new (a[th2fSize]) TH2F(name, title, bins1, xlow1, xhigh1, bins2, xlow2, xhigh2);
//  std::cout << "histoTH2F    " << name << std::endl;
  th2fSize++;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//Add new TH2F histogram
//
/////////////////////////////////////////////////////////////////////////////////////
void HistoManager::addTH2F(const char *name, const char *title, int binsx, const double *xbins , int binsy, const double *ybins) {
  TClonesArray &a = *m_TH2FArray;
  a[th2fSize] = new (a[th2fSize]) TH2F(name, title, binsx, xbins, binsy, ybins);
//  std::cout << "histoTH2F    " << name << std::endl;
  th2fSize++;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//Add new TProfile histogram
//
/////////////////////////////////////////////////////////////////////////////////////
void HistoManager::addTProf(const char *name, const char *title, int bins, float xlow1, float xhigh1,float xlow2, float xhigh2) {
  TClonesArray &a = *m_TProfArray;
  a[tprofSize] = new (a[tprofSize]) TProfile(name, title, bins, xlow1, xhigh1, xlow2, xhigh2);
//  std::cout << "histoTProfile    " << name << std::endl;
  tprofSize++;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//Create all defined  histograms 
//
/////////////////////////////////////////////////////////////////////////////////////
void HistoManager::CreateHistTables() {

   //  Create tables which host histograms...

  //  TH1::SetDefaultSumw2(kTrue);
  //  TH1::SetDefaultSumw2(kFalse);


  //-------------------------------------//
  m_TH1FArray = new TClonesArray("TH1F",5);

  //-------------------------------------//
  m_TH2FArray = new TClonesArray("TH2F",8);

  //-------------------------------------//
  m_TH3FArray = new TClonesArray("TH3F",11);

  //-------------------------------------//
  m_TProfArray = new TClonesArray("TProfile",7);

}

/////////////////////////////////////////////////////////////////////////////////////
TProfile *HistoManager::GetHistoTProf(int x){ // returns the TProfile specified

  TProfile *histo_prof;

  Text_t name[500];

  sprintf(name,"hist%.2d",x);

  histo_prof = (TProfile*) m_TProfArray->FindObject(name);
  if (histo_prof == NULL)
    std::cout << "NULL histogram hist" <<x <<  std::endl;

  return (TProfile*) m_TProfArray->FindObject(name);
}

/////////////////////////////////////////////////////////////////////////////////////
TH1F *HistoManager::GetHistoTH1F(int x){// returns the histogram specified
  
  TH1F *histo_f;
  
  char name[500];
  sprintf(name,"hist%.2d",x);
  
  histo_f =  (TH1F*)m_TH1FArray->FindObject(name);
  if (histo_f == NULL)
    std::cout << "NULL histogram hist" <<x <<  std::endl;

  return histo_f;
}

/////////////////////////////////////////////////////////////////////////////////////
TH2F *HistoManager::GetHistoTH2F(int x){// returns the histogram specified
  
  TH2F *histo_f;
  
  char name[500];
  sprintf(name,"histomap%.2d",x);
  
  histo_f =  (TH2F*)m_TH2FArray->FindObject(name);
  if (histo_f == NULL)
    std::cout << "NULL histomap " << x << std::endl;

  return histo_f;
}

/////////////////////////////////////////////////////////////////////////////////////
TH3F *HistoManager::GetHistoTH3F(int x){// returns the histogram specified
  
  TH3F *histo_f;
  
  char name[500];
  sprintf(name,"histomap%.2d",x);
  
  histo_f =  (TH3F*)m_TH3FArray->FindObject(name);
  if (histo_f == NULL)
    std::cout << "NULL histomap " << x << std::endl;

  return histo_f;
}
  
/////////////////////////////////////////////////////////////////////////////////////
void HistoManager::StoreHistos(){ // saves all the histograms in the TClonesArray 
  
  m_TH1FArray->Write();
  m_TH2FArray->Write();
  m_TH3FArray->Write();
  m_TProfArray->Write();

}


