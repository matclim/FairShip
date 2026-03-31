#ifndef LASTBITMUONSHIELDHIT_H
#define LASTBITMUONSHIELDHIT_H 1
#include "FairVolume.h"
#include "ShipHit.h"
#include "lastBitMuonShieldPoint.h"
#include "TObject.h"
#include "TGeoShape.h"
#include "TGeoPhysicalNode.h"


class lastBitMuonShieldHit : public ShipHit {
public:
  
  /** Default constructor **/
  lastBitMuonShieldHit();

  /** Constructor from lastBitMuonShieldHit
   *@param detID    Detector ID
   *@param t_1, t_2      TDC on both sides
   *@param flag      True/False, in case of pile up
   **/
  lastBitMuonShieldHit(lastBitMuonShieldPoint* p, Double_t t0);

  /** Destructor **/
  virtual ~lastBitMuonShieldHit();
  
  /** Accessors **/
  Double_t GetX();
  Double_t GetY();
  Double_t GetZ();
  TVector3 GetXYZ();
  TGeoNode* GetNode();
  std::vector<double> GetTime(Double_t x);
  std::vector<double> GetTime();
  std::vector<double> GetMeasurements();
  /** Modifier **/
  void SetTDC(Float_t val1, Float_t val2){t_1=val1;t_2=val2;}
  
  /** Output to screen **/
  virtual void Print() const;

  void Dist(Float_t x, Float_t& lpos, Float_t& lneg);
  Double_t Resol(Double_t x);
  void setInvalid() {flag = false;}
  void setIsValid() {flag = true;}
  bool isValid() const {return flag;}
private:
  lastBitMuonShieldHit(const lastBitMuonShieldHit& point);
  lastBitMuonShieldHit operator=(const lastBitMuonShieldHit& point);
  Double_t v_drift = 15.; // cm/ns
  Double_t par[4] = { 0.0272814, 109.303, 0, 0.0539487 };
  
  Float_t flag;     ///< flag
  Float_t t_1,t_2;  ///< TDC on both sides
  
  ClassDef(lastBitMuonShieldHit,1);

};

#endif
