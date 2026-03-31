#ifndef UPSTREAMTAGGER_H
#define UPSTREAMTAGGER_H

#include "FairDetector.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "ShipUnit.h"

class UpstreamTaggerPoint;
class FairVolume;
class TClonesArray;

using ShipUnit::m;
using ShipUnit::cm;

class UpstreamTagger: public FairDetector
{

  public:

    /**      Name :  Detector Name
     *       Active: kTRUE for active detectors (ProcessHits() will be called)
     *               kFALSE for inactive detectors
    */
    UpstreamTagger(const char* Name, Bool_t Active);

    /** default constructor */
    UpstreamTagger();

    /** destructor */
    virtual ~UpstreamTagger();

    /** Initialization of the detector is done here */
    virtual void   Initialize();

    /**   this method is called for each step during simulation
     *    (see FairMCApplication::Stepping())
    */
    virtual Bool_t ProcessHits( FairVolume* v=0);

    /**       Registers the produced collections in FAIRRootManager. */
    virtual void Register();

    /** Gets the produced collections */
    virtual TClonesArray* GetCollection(Int_t iColl) const;

    /** has to be called after each event to reset the containers */
    virtual void Reset();

  /** Sets detector dimensions */
  void SetXdim(Double_t x) {det_xdim = x;}
  void SetYdim(Double_t y) {det_ydim = y;}
  void SetZdim(Double_t z) {det_zdim = z;}
  
    /** Sets detector position along z */
    void SetZposition(Double_t z) {det_zPos = z;}

    /**  Create the detector geometry */
    void ConstructGeometry();

    /**      This method is an example of how to add your own point
     *       of type TimeRpcPoint to the clones array
    */
    UpstreamTaggerPoint* AddHit(Int_t trackID, Int_t detID,
			 TVector3 pos, TVector3 mom,
			 Double_t time, Double_t length,
			 Double_t eLoss, Int_t pdgCode,TVector3 Lpos, TVector3 Lmom);

    virtual void   EndOfEvent();
    virtual void   FinishPrimary() {;}
    virtual void   FinishRun() {;}
    virtual void   BeginPrimary() {;}
    virtual void   PostTrack() {;}
    virtual void   PreTrack() {;}
    virtual void   BeginEvent() {;}

    /** Track information to be stored until the track leaves the active volume.*/
    Int_t          fTrackID;            //!  track index
    Int_t          fVolumeID;           //!  volume id
    TLorentzVector fPos;                //!  position at entrance
    TLorentzVector fMom;                //!  momentum at entrance
    Double_t       fTime;               //!  time
    Double_t       fLength;             //!  length
    Double_t       fELoss;              //!  energy loss

  /** Detector parameters.*/
  Double_t det_xdim; //!  x-dim of veto station 
  Double_t det_ydim; //!  y-dim of veto station
  Double_t det_zdim; //!  z-dim of veto station

  Double_t det_zPos;     //!  z-position of veto station

  private:

    TGeoVolume* UpstreamTagger_fulldet; // Timing_detector_1 object

    /** container for data points */
    TClonesArray* fUpstreamTaggerPointCollection;

    UpstreamTagger(const UpstreamTagger&);
    UpstreamTagger& operator=(const UpstreamTagger&);
    Int_t InitMedium(const char* name);


    ClassDef(UpstreamTagger,1)
};

#endif //UPSTREAMTAGGER_H
