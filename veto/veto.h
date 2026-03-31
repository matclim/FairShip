#ifndef VETO_VETO_H_
#define VETO_VETO_H_

#include "FairDetector.h"
#include "TGeoVolume.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <map>

class vetoPoint;
class FairVolume;
class TClonesArray;

/*
 *
 * @class Veto
 * @brief Class representing the Veto detector.
 *
 * The Veto class simulates the geometry and behavior of the veto detector
 * within the FairShip framework. It includes methods for initializing
 * the detector geometry, processing hits, and handling input parameters.
 *
*/

class veto : public FairDetector {
public:
  /*
   *      Name: Veto
   *      Active: kTRUE for active detectors (ProcessHits() will be called)
   *              kFALSE for inactive detectors
   */
  veto();
  /*
   *       destructor     
   */
  virtual ~veto();
  /*
   *      Initialization of the detector is done here    
   */
  virtual void Initialize();
  /*
   *       this method is called for each step during simulation
   *       (see FairMCApplication::Stepping())
   */
  virtual Bool_t ProcessHits(FairVolume* v = 0);
  /*
   *       Registers the produced collections in FAIRRootManager.     
   */
  virtual void Register();
  /*
   * Gets the produced collections 
   */
  virtual TClonesArray* GetCollection(Int_t iColl) const;
  /*
   *      has to be called after each event to reset the containers      
   */
  virtual void Reset();
  
  void SetFastMuon() { fFastMuon = true; }       // kill all tracks except of muons
  void SetFollowMuon() { fFollowMuon = true; }   // make muon shield active to follow muons
  void SetVesselDimensions(Float_t x1, Float_t x2, Float_t y1, Float_t y2, Float_t vesselstart) {
    VetoStartInnerX = x1;
    VetoEndInnerX = x2;
    VetoStartInnerY = y1;
    VetoEndInnerY = y2;
    zStartDecayVol = vesselstart;
  }
  /*
   *      Create the detector geometry        
   */
  void ConstructGeometry();
  void SetVesselStructure(Double_t a,
			  TString  b,
			  TString  c) {
    f_VetoThickness = a;
    vetoMed_name = b;
    decayVolumeMed_name = c;
  }

  /*
   *      This method is an example of how to add your own point
   *       of type vetoPoint to the clones array
   */
  vetoPoint* AddHit(Int_t trackID,
		    Int_t detID,
		    TVector3 pos,
		    TVector3 mom,
		    Double_t time,
		    Double_t length,
		    Double_t eLoss,
		    Int_t pdgcode,
		    TVector3 Lpos,
		    TVector3 Lmom);

  /*
   * The following methods can be implemented if you need to make
   *  any optional action in your detector during the transport.
   */
  
  virtual void CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset) { ; }
  virtual void SetSpecialPhysicsCuts() { ; }
  virtual void EndOfEvent();
  virtual void FinishPrimary() { ; }
  virtual void FinishRun() { ; }
  virtual void BeginPrimary() { ; }
  virtual void PostTrack() { ; }
  virtual void PreTrack();
  virtual void BeginEvent() { ; }
  
private:
  /*
   * Track information to be stored until the track leaves the active volume.
   */
  //!  track index
  Int_t fTrackID;
  //!  volume id
  Int_t fVolumeID;
  //!  position at entrance
  TLorentzVector fPos;
  //!  momentum at entrance
  TLorentzVector fMom;
  //!  time
  Float_t fTime;
  //!  length
  Float_t fLength;
  //!  energy loss
  Float_t fELoss;
  
  Bool_t fFastMuon, fFollowMuon;

  //! thickness of veto detector
  Double_t f_VetoThickness;
  //! medium of veto counter (default: vacuum) 
  TString vetoMed_name;
  //! medium of decay volume (default: helium)
  TString decayVolumeMed_name;
  
  TGeoMedium* vetoMed;
  TGeoMedium* decayVolumeMed;
  
  //! Width of the Vessel along X at the start
  Float_t VetoStartInnerX;
  //! Length of the Vessel along Y at the start
  Float_t VetoStartInnerY;
  //! Width of the Vessel along X at the end
  Float_t VetoEndInnerX;
  //! Length of the Vessel along Y at the end
  Float_t VetoEndInnerY;
  
  //! z Position of the Decay Volume start in the global coordinate system
  Float_t zStartDecayVol;
  
  /** container for data points */
  TClonesArray* fvetoPointCollection;
  
  veto(const veto&);
  veto& operator=(const veto&);
  Int_t InitMedium(const char* name);
  /*
   * Adds a solid Trapezoid of thickness (along z) wz with start cross-section dimensions of wX_start * wY_start
   * and end cross-section dimensions of wX_end *wY_end
   */
  TGeoVolume* GeoTrapezoid(TString xname,
			   Double_t wz,
			   Double_t wX_start,
			   Double_t wX_end,
			   Double_t wY_start,
			   Double_t wY_end,
			   Int_t color,
			   TGeoMedium* material,
			   Bool_t sens);

  /*
   * Adds a Hollow Trapezoid of thickness (along z) wz with start cross-section dimensions of wX_start * wY_start
   * and end cross-section dimensions of wX_end *wY_end. The trapezoid is hollowed out to have a thickness of "thick" cm.
   */
  TGeoVolume* GeoTrapezoidHollow(TString xname,
				 Double_t thick,
				 Double_t wz,
				 Double_t wX_start,
				 Double_t wX_end,
				 Double_t wY_start,
				 Double_t wY_end,
				 Int_t color,
				 TGeoMedium* material,
				 Bool_t sens);
  
  /*
   * Adds a custom block of OuterWall+InnerWall ribs for a given distance along z.
   * Ensures consistency in implementation throught the z.
   */
  void AddBlock(TGeoVolumeAssembly* tInnerWall,
		TGeoVolumeAssembly* tDecayVacuum,
		int blockNr,
		double z1,
		double z2,
		double Zshift,
		double wallThick);
  
  //! slope along the width (x)
  double wx(double z);
  //! slope along the length (y)
  double wy(double z);
  
  TGeoVolume* MakeSegments();
  
  ClassDef(veto, 1)
};

#endif   // VETO_VETO_H_
