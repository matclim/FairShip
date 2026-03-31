#include "lastBitMuonShield.h"
#include "lastBitMuonShieldPoint.h"

#include "FairVolume.h"
#include "FairGeoVolume.h"
#include "FairGeoNode.h"
#include "FairRootManager.h"
#include "FairGeoLoader.h"
#include "FairGeoInterface.h"
#include "FairGeoMedia.h"
#include "FairGeoBuilder.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "ShipDetectorList.h"
#include "ShipStack.h"

#include "TClonesArray.h"
#include "TVirtualMC.h"
#include "TGeoManager.h"
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoTube.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TMath.h"
#include "TParticle.h"
#include "TVector3.h"

#include <iostream>
#include <sstream>
using std::cout;
using std::endl;


lastBitMuonShield::lastBitMuonShield() : FairDetector("lastBitMuonShield", kTRUE, klastBitMuonShield),
					 fTrackID(-1),
					 fVolumeID(-1),
					 fPos(),
					 fMom(),
					 fTime(-1.),
					 fLength(-1.),
					 fELoss(-1),
					 //
					 det_xdim(0),
					 det_ydim(0),
					 det_zdim(0),
					 det_zPos(0),
					 //
					 flastBitMuonShieldPointCollection(new TClonesArray("lastBitMuonShieldPoint")) {
}



lastBitMuonShield::lastBitMuonShield(const char* name, Bool_t active) : FairDetector(name, active, klastBitMuonShield),
									fTrackID(-1),
									fVolumeID(-1),
									fPos(),
									fMom(),
									fTime(-1.),
									fLength(-1.),
									fELoss(-1),
									//
									det_xdim(0),
									det_ydim(0),
									det_zdim(0),
									det_zPos(0),
									//
									flastBitMuonShieldPointCollection(new TClonesArray("lastBitMuonShieldPoint")) {
}

void lastBitMuonShield::Initialize() {
  FairDetector::Initialize();
}


lastBitMuonShield::~lastBitMuonShield() {
  if (flastBitMuonShieldPointCollection) {
    flastBitMuonShieldPointCollection->Delete();
    delete flastBitMuonShieldPointCollection;
  }
}

Int_t lastBitMuonShield::InitMedium(const char* name) {
   static FairGeoLoader *geoLoad=FairGeoLoader::Instance();
   static FairGeoInterface *geoFace=geoLoad->getGeoInterface();
   static FairGeoMedia *media=geoFace->getMedia();
   static FairGeoBuilder *geoBuild=geoLoad->getGeoBuilder();

   FairGeoMedium *ShipMedium=media->getMedium(name);

   if (!ShipMedium) {
     Fatal("InitMedium","Material %s not defined in media file.", name);
     return -1111;
   }
   
   TGeoMedium* medium=gGeoManager->GetMedium(name);
   if (medium!=NULL) return ShipMedium->getMediumIndex();

   return geoBuild->createMedium(ShipMedium);

  return 0;
}

Bool_t  lastBitMuonShield::ProcessHits(FairVolume* vol) {
  /** This method is called from the MC stepping */
  //Set parameters at entrance of volume. Reset ELoss.
  if ( gMC->IsTrackEntering() ) {
    fELoss  = 0.;
    fTime   = gMC->TrackTime() * 1.0e09;
    fLength = gMC->TrackLength();
    gMC->TrackPosition(fPos);
    gMC->TrackMomentum(fMom);
  }

  // Sum energy loss for all steps in the active volume
  fELoss += gMC->Edep();

  // Create vetoPoint at exit of active volume
  if ( gMC->IsTrackExiting()    ||
       gMC->IsTrackStop()       ||
       gMC->IsTrackDisappeared()   ) {
    if (fELoss == 0. ) { return kFALSE; }
    
    fTrackID  = gMC->GetStack()->GetCurrentTrackNumber();
    
    Int_t uniqueId;
    gMC->CurrentVolID(uniqueId);
    if (uniqueId>1000000) { //Solid scintillator case {
      Int_t vcpy;
      gMC->CurrentVolOffID(1, vcpy);
      if (vcpy==5) uniqueId+=4; //Copy of half
    }
    
    TParticle* p = gMC->GetStack()->GetCurrentTrack();
    Int_t pdgCode = p->GetPdgCode();
    TLorentzVector Pos;
    gMC->TrackPosition(Pos);
    TLorentzVector Mom;
    gMC->TrackMomentum(Mom);
    Double_t xmean = (fPos.X()+Pos.X())/2. ;
    Double_t ymean = (fPos.Y()+Pos.Y())/2. ;
    Double_t zmean = (fPos.Z()+Pos.Z())/2. ;
    
    //cout << "MELQUI DE PRUSIA" << " :(" << xmean << ", " << ymean << ", " << zmean << "): " << gMC->CurrentVolName() << endl;
    
    AddHit(fTrackID, uniqueId, TVector3(xmean, ymean,  zmean),
           TVector3(fMom.Px(), fMom.Py(), fMom.Pz()), fTime, fLength,
           fELoss,pdgCode,TVector3(Pos.X(), Pos.Y(), Pos.Z()),
	   TVector3(Mom.Px(), Mom.Py(), Mom.Pz()) );
    
    // Increment number of veto det points in TParticle
    ShipStack* stack = (ShipStack*) gMC->GetStack();
    stack->AddPoint(klastBitMuonShield);
  }
  
  return kTRUE;
}

void lastBitMuonShield::EndOfEvent() {
  flastBitMuonShieldPointCollection->Clear();
}



void lastBitMuonShield::Register(){

  /** This will create a branch in the output tree called
      lastBitMuonShieldPoint, setting the last parameter to kFALSE means:
      this collection will not be written to the file, it will exist
      only during the simulation.
  */

  FairRootManager::Instance()->Register("lastBitMuonShieldPoint", "lastBitMuonShield",
                                        flastBitMuonShieldPointCollection, kTRUE);
}



TClonesArray* lastBitMuonShield::GetCollection(Int_t iColl) const {
  if (iColl == 0) { return flastBitMuonShieldPointCollection; }
  else { return NULL; }
}



void lastBitMuonShield::Reset()
{
  flastBitMuonShieldPointCollection->Clear();
}



void lastBitMuonShield::ConstructGeometry() {
  TGeoVolume *top = gGeoManager->GetTopVolume();
  
  InitMedium("iron");
  TGeoMedium *iron =gGeoManager->GetMedium("iron");
  
  //////////////////////////////////////////////////////
  auto *plate = new TGeoBBox("lastBitMuonShield",
			     det_xdim/2.0,
			     det_ydim/2.0,
			     det_zdim/2.0);
  auto *plate_geo_vol = new TGeoVolume("lastBitMuonShield_geo_vol", plate, iron);
  plate_geo_vol -> SetLineColor(kBlue);
  AddSensitiveVolume(plate_geo_vol);
  top->AddNode(plate_geo_vol, 1, new TGeoTranslation(0, 0, det_zPos));

  ///////////////////////////////////////////////////////

  return;
}



lastBitMuonShieldPoint* lastBitMuonShield::AddHit(Int_t trackID, Int_t detID,
						  TVector3 pos, TVector3 mom,
						  Double_t time, Double_t length,
						  Double_t eLoss, Int_t pdgCode,TVector3 Lpos, TVector3 Lmom) {
  TClonesArray& clref = *flastBitMuonShieldPointCollection;
  Int_t size = clref.GetEntriesFast();
  // cout << "veto hit called "<< pos.z()<<endl;
  return new(clref[size]) lastBitMuonShieldPoint(trackID, detID, pos, mom,
		         time, length, eLoss, pdgCode,Lpos,Lmom);
}
