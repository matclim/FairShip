// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP
// Collaboration

#ifndef UPSTREAMTAGGER_UPSTREAMTAGGERHIT_H_
#define UPSTREAMTAGGER_UPSTREAMTAGGERHIT_H_

#include "DetectorHit.h"
#include "TVector3.h"

class UpstreamTaggerPoint;

/**
 * @brief Hit class for UpstreamTagger scoring plane
 *
 * Simple hit class for UBT scoring plane detector.
 * Stores smeared position and time from MC truth.
 * Does not store MC truth information directly.
 */
class UpstreamTaggerHit : public SHiP::DetectorHit {
 public:
  /** Default constructor **/
  UpstreamTaggerHit();

  /** Constructor from UpstreamTaggerPoint
   * @param p     MC point
   * @param t0    Event time offset
   * @param pos_res Position resolution (cm)
   * @param time_res Time resolution (ns)
   **/
  UpstreamTaggerHit(UpstreamTaggerPoint* p, Double_t t0, Double_t pos_res,
                    Double_t time_res);

  /** Destructor **/
  ~UpstreamTaggerHit() override = default;

  /** Copy constructor **/
  UpstreamTaggerHit(const UpstreamTaggerHit& hit) = default;
  UpstreamTaggerHit& operator=(const UpstreamTaggerHit& hit) = default;

<<<<<<< HEAD
  /** Position accessors **/
  Double_t GetX() const { return fX; }
  Double_t GetY() const { return fY; }
  Double_t GetZ() const { return fZ; }
  TVector3 GetXYZ() const { return TVector3(fX, fY, fZ); }

  /** Time accessor **/
  Double_t GetTime() const { return fTime; }

  /** Output to screen **/
  using SHiP::DetectorHit::Print;
  void Print() const;

 private:
  Double_t fX;     ///< Smeared x position (cm)
  Double_t fY;     ///< Smeared y position (cm)
  Double_t fZ;     ///< Smeared z position (cm)
  Double_t fTime;  ///< Smeared time (ns)
=======
    TGeoNode* GetNode(Double_t &hit_final, Int_t &mod);
    std::vector<double> GetTime(Double_t x);
    std::vector<double> GetTime();
    std::vector<double> GetMeasurements();
    /** Modifier **/
    void SetPoint(Double_t p1, Double_t p2, Double_t p3){point_final[0]=p1;point_final[1]=p2;point_final[2]=p3;}
    /** Output to screen **/
    virtual void Print() const;

    void setInvalid() {flag = false;}
    void setIsValid() {flag = true;}

    //Rpc time is invalid if isValid returns False
    bool isValid() const {return flag;}
  private:
    UpstreamTaggerHit(const UpstreamTaggerHit& point);
    UpstreamTaggerHit operator=(const UpstreamTaggerHit& point);

    UpstreamTagger* c0;
    Double_t point_final[3];
    const Double_t * mom[3];

    Float_t flag;     ///< flag
    Float_t time;
    Double_t X, Y, Z;
  
    ClassDef(UpstreamTaggerHit,1);
>>>>>>> 773ac29b8 (LastBitMuonShield, dummyUBT and cave set to vacuum)

  ClassDef(UpstreamTaggerHit, 2);
};

#endif  // UPSTREAMTAGGER_UPSTREAMTAGGERHIT_H_
