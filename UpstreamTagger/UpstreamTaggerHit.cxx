#include "UpstreamTaggerHit.h"
#include "UpstreamTagger.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRandom1.h"
#include "TRandom3.h"
#include "TGeoManager.h"
#include "TGeoBBox.h"
#include "TGeoNode.h"
#include "TRandom3.h"

#include <iostream>
#include <cmath>
#include <stdlib.h>     /* srand, rand */
#include <cstdlib>
#include <ctime>
#include <time.h>       /* time */

using std::cout;
using std::endl;

// -----   Default constructor   --------------
UpstreamTaggerHit::UpstreamTaggerHit()
  : ShipHit()
{
 flag = true;
}


// -----   constructor from TimeDetPoint from TimeDetHit-------------------------------
UpstreamTaggerHit::UpstreamTaggerHit(UpstreamTaggerPoint* p, UpstreamTagger* c, Double_t t0)
  : ShipHit()
{

     fDetectorID = p->GetDetectorID();
     c0 = c;

     Double_t a1, b1, c1;
     a1 = p->GetX(); b1 = p->GetY(); c1 = p->GetZ();
     const Double_t * point1[3];
     point1[0] = &a1; point1[1] = &b1; point1[2] = &c1;
     point_final[0] = *point1[0];point_final[1] = *point1[1];point_final[2] = *point1[2];

     Double_t a2, b2, c2;
     a2 = p->GetPx(); b2 = p->GetPy(); c2 = p->GetPz();
     mom[0] = &a2; mom[1] = &b2; mom[2] = &c2;

     time = p->GetTime();
}


// -----   Destructor   -------------------------
UpstreamTaggerHit::~UpstreamTaggerHit() { }

// ---- return time information for a given track extrapolation
std::vector<double>  UpstreamTaggerHit::GetTime(Double_t x){
     std::vector<double> retVal;
     retVal.push_back(time);
     return retVal;
}
// ---- return mean time information
std::vector<double>  UpstreamTaggerHit::GetTime(){
     TGeoBBox* shape =  (TGeoBBox*)gGeoManager->GetVolume("UpstreamTagger")->GetShape();
   
     std::vector<double> retVal;
     retVal.push_back(time);
     return retVal;
}

// ----------------------------------------------
TVector3 UpstreamTaggerHit::GetXYZ()
{
    Double_t hit_final;
    Int_t mod;
    TGeoNavigator* nav = gGeoManager->GetCurrentNavigator();
    TGeoNode* node = GetNode(hit_final, mod);
    auto shape =  dynamic_cast<TGeoBBox*>(node->GetVolume()->GetShape());
    Double_t origin[3] = {shape->GetOrigin()[0],shape->GetOrigin()[1],shape->GetOrigin()[2]};
    Double_t master[3] = {0,0,0};
    nav->LocalToMaster(origin,master);
    TVector3 pos = TVector3(master[0],master[1],master[2]);
    return pos;
}


Double_t UpstreamTaggerHit::GetX()
{ TVector3 pos = GetXYZ();
  return pos.X();
}


Double_t UpstreamTaggerHit::GetY()
{ TVector3 pos = GetXYZ();
  return pos.Y();
}


Double_t UpstreamTaggerHit::GetZ()
{ TVector3 pos = GetXYZ();
  return pos.Z();
}


TGeoNode* UpstreamTaggerHit::GetNode(Double_t &hit_final, Int_t &mod) {
  TGeoNavigator* nav = gGeoManager->GetCurrentNavigator();
  mod = 0;

  Double_t XHit1 = point_final[0];
  Double_t YHit1 = point_final[1];
  Double_t ZHit1 = point_final[2];
  
  const Double_t PxHit1 = *mom[0];
  const Double_t PyHit1 = *mom[1];
  const Double_t PzHit1 = *mom[2];
  
  Z = ZHit1;
  X = XHit1 + (Z*(PxHit1/TMath::Sqrt(PxHit1*PxHit1 + PyHit1*PyHit1 + PzHit1*PzHit1)));
  Y = YHit1 + (Z*(PyHit1/TMath::Sqrt(PxHit1*PxHit1 + PyHit1*PyHit1 + PzHit1*PzHit1)));
  
  hit_final = X;
  hit_final = Y;

  c0->det_zPos;
  
  return nav->GetCurrentNode();
}

// -----   Public method Print   -----------------------
void UpstreamTaggerHit::Print() const
{
  cout << "-I- UpstreamTaggerHit: UpstreamTagger hit " << " in detector " << fDetectorID << endl;
  cout << "MIMI (det_zPos):" << det_zPos << "HIT (X): " << X << "HIT (Y): " << Y << "HIT (Z): " << Z << endl;
}


// -----------------------------------------------------
