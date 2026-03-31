
// 09/07/2024
// SBT software contact: anupama.reghunath@cern.ch
/*
 *
 * @file veto.cxx
 * @brief Implementation of the Veto detector class.
 *
 * This file contains the definitions for the Veto class used in the
 * FairShip Software. The class is responsible for simulating the Decay Vessel(helium) + integrated SBT geometry and
 * interactions.
 */

#include "veto.h"

#include "FairGeoBuilder.h"
#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoMedia.h"
#include "FairGeoNode.h"
#include "FairGeoVolume.h"
#include "FairLogger.h"   /// for FairLogger, MESSAGE_ORIGIN
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "FairVolume.h"
#include "ShipDetectorList.h"
#include "ShipStack.h"
#include "TClonesArray.h"
#include "TGeoArb8.h"
#include "TGeoBBox.h"
#include "TGeoBoolNode.h"
#include "TGeoCompositeShape.h"
#include "TGeoCone.h"
#include "TGeoEltu.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoShapeAssembly.h"
#include "TGeoSphere.h"
#include "TGeoTube.h"
#include "TMath.h"
#include "TParticle.h"
#include "TVirtualMC.h"
#include "vetoPoint.h"

#include <iostream>
#include <math.h>
#include <vector>

Double_t cm = 1;          // cm
Double_t m = 100 * cm;    //  m
Double_t mm = 0.1 * cm;   //  mm

/**
 * @brief Constructor for the Veto class.
 *
 * Initializes the veto detector with the given parameters.
 *
 */
veto::veto()
    : FairDetector("Veto", kTRUE, kVETO)
    , fTrackID(-1)
    , fVolumeID(-1)
    , fPos()
    , fMom()
    , fTime(-1.)
    , fLength(-1.)
    , fELoss(-1)
    , fvetoPointCollection(new TClonesArray("vetoPoint"))
    , fFastMuon(kFALSE)
    , fFollowMuon(kFALSE) {}

/*
 *
 * @brief Destructor for the Veto class.
 *
 * Cleans up any resources used by the Veto detector.
 */

veto::~veto() {
  if (fvetoPointCollection) {
    fvetoPointCollection->Delete();
    delete fvetoPointCollection;
  }
}

void veto::Initialize() {
  FairDetector::Initialize();
}

TGeoVolume* veto::GeoTrapezoid(TString xname,
                               Double_t z_thick,
                               Double_t x_thick_start,
                               Double_t x_thick_end,
                               Double_t y_thick_start,
                               Double_t y_thick_end,
                               Int_t color,
                               TGeoMedium* material,
                               Bool_t sens = kFALSE) {
  Double_t dz = z_thick / 2;
  Double_t dx1 = x_thick_start / 2 - 1.E-6;
  Double_t dx2 = x_thick_end / 2 - 1.E-6;
  Double_t dy1 = y_thick_start / 2 - 1.E-6;
  Double_t dy2 = y_thick_end / 2 - 1.E-6;
  TGeoArb8* T1 = new TGeoArb8("T1" + xname, dz + 1.E-6);
  T1->SetVertex(0, -dx1, -dy1);
  T1->SetVertex(1, -dx1, dy1);
  T1->SetVertex(2, dx1, dy1);
  T1->SetVertex(3, dx1, -dy1);
  T1->SetVertex(4, -dx2, -dy2);
  T1->SetVertex(5, -dx2, dy2);
  T1->SetVertex(6, dx2, dy2);
  T1->SetVertex(7, dx2, -dy2);
  
  TGeoVolume* T = new TGeoVolume(xname, T1, material);
  T->SetLineColor(color);
  if (sens) {
    AddSensitiveVolume(T);
  }   // and make the volumes sensitive..
  return T;
}

TGeoVolume* veto::GeoTrapezoidHollow(TString xname,
                                     Double_t wallthick,
                                     Double_t z_thick,
                                     Double_t x_thick_start,
                                     Double_t x_thick_end,
                                     Double_t y_thick_start,
                                     Double_t y_thick_end,
                                     Int_t color,
                                     TGeoMedium* material,
                                     Bool_t sens = kFALSE) {
  
  Double_t dx_start = x_thick_start / 2;
  Double_t dy_start = y_thick_start / 2;
  Double_t dx_end = x_thick_end / 2;
  Double_t dy_end = y_thick_end / 2;
  Double_t dz = z_thick / 2;
  
  TString nm = xname.ReplaceAll("-", "");   // otherwise it will try to subtract "-" in TGeoComposteShape
  Double_t dx1 = dx_start + wallthick;
  Double_t dx2 = dx_end + wallthick;
  Double_t dy1 = dy_start + wallthick;
  Double_t dy2 = dy_end + wallthick;
  
  TGeoArb8* T2 = new TGeoArb8("T2" + nm, dz);
  T2->SetVertex(0, -dx1, -dy1);
  T2->SetVertex(1, -dx1, dy1);
  T2->SetVertex(2, dx1, dy1);
  T2->SetVertex(3, dx1, -dy1);
  T2->SetVertex(4, -dx2, -dy2);
  T2->SetVertex(5, -dx2, dy2);
  T2->SetVertex(6, dx2, dy2);
  T2->SetVertex(7, dx2, -dy2);
  
  Double_t tdx1 = dx_start;
  Double_t tdx2 = dx_end;
  Double_t tdy1 = dy_start;
  Double_t tdy2 = dy_end;
  TGeoArb8* T1 = new TGeoArb8("T1" + nm, dz + 1.E-6);
  T1->SetVertex(0, -tdx1, -tdy1);
  T1->SetVertex(1, -tdx1, tdy1);
  T1->SetVertex(2, tdx1, tdy1);
  T1->SetVertex(3, tdx1, -tdy1);
  T1->SetVertex(4, -tdx2, -tdy2);
  T1->SetVertex(5, -tdx2, tdy2);
  T1->SetVertex(6, tdx2, tdy2);
  T1->SetVertex(7, tdx2, -tdy2);
  
  TGeoCompositeShape* T321 = new TGeoCompositeShape("T3" + nm, "T2" + nm + "-T1" + nm);
  TGeoVolume* T = new TGeoVolume(xname, T321, material);
  T->SetLineColor(color);
  // and make the volumes sensitive..
  if (sens) {
    AddSensitiveVolume(T);
  }
  return T;
}

Double_t veto::wx(Double_t z) {
  // calculate x thickness at z
  Double_t wx1 = VetoStartInnerX;
  Double_t wx2 = VetoEndInnerX;
  Double_t z1 = 0 * m;
  Double_t z2 = 50 * m;
  return (wx1 + (z - z1) * (wx2 - wx1) / (z2 - z1));
}

Double_t veto::wy(Double_t z) {
  // calculate y thickness at z
  Double_t wy1 = VetoStartInnerY;
  Double_t wy2 = VetoEndInnerY;
  Double_t z1 = 0 * m;
  Double_t z2 = 50 * m;
  return (wy1 + (z - z1) * (wy2 - wy1) / (z2 - z1));
}

void veto::AddBlock(TGeoVolumeAssembly* tInnerWall,
                    TGeoVolumeAssembly* tDecayVacuum,
		    int blockNr,
                    Double_t z1,
                    Double_t z2,
                    Double_t Zshift,
                    Double_t wallThick) {
  TString blockName = "block";
  blockName += blockNr;
  
  int ribColor = 15;
  Double_t wz = (z2 - z1);

  // inner wall
  TString nameInnerWall = (TString) tInnerWall -> GetName() + "_" + blockName;

  TGeoVolume* TIW = GeoTrapezoidHollow(nameInnerWall,
				       wallThick,
				       wz,
				       wx(z1),
				       wx(z2),
				       wy(z1),
				       wy(z2),
				       ribColor,
				       vetoMed);
  
  tInnerWall->AddNode(TIW, 0, new TGeoTranslation(0, 0, Zshift));

  // decay vacuum
  TString nameDecayVacuum = (TString) tDecayVacuum -> GetName() + "_" + blockName;
    
  TGeoVolume* TDV = GeoTrapezoid(nameDecayVacuum,
				 wz,
				 wx(z1),
				 wx(z2),
				 wy(z1),
				 wy(z2),
				 1,
				 decayVolumeMed);

  TDV->SetVisibility(kFALSE);
  tDecayVacuum->AddNode(TDV, 0, new TGeoTranslation(0, 0, Zshift));
}

TGeoVolume* veto::MakeSegments() {
  TGeoVolumeAssembly* tTankVol = new TGeoVolumeAssembly("T2");

  TString nameInnerWall = "VetoInnerWall";
  TGeoVolumeAssembly* tInnerWall = new TGeoVolumeAssembly(nameInnerWall);
  
  TString nameDecayVacuum = "DecayVacuum";
  TGeoVolumeAssembly* tDecayVacuum = new TGeoVolumeAssembly(nameDecayVacuum);

  Double_t wallThick = f_VetoThickness;
  
  Double_t z1 = 0 * m;
  Double_t z2 = 50.0 * m;
  Double_t wz = (z2 - z1);
    
  Double_t Zshift = wz / 2;   // calibration of Z position
  
  AddBlock(tInnerWall,
	   tDecayVacuum,
	   1,
	   z1,
	   z2,
	   Zshift,
	   wallThick);  
  
  tTankVol->AddNode(tInnerWall, 0, new TGeoTranslation(0, 0, 0));
  tTankVol->AddNode(tDecayVacuum, 0, new TGeoTranslation(0, 0, 0));
    
  return tTankVol;
}

// -----   Private method InitMedium
Int_t veto::InitMedium(const char* name) {
  static FairGeoLoader* geoLoad = FairGeoLoader::Instance();
  static FairGeoInterface* geoFace = geoLoad->getGeoInterface();
  static FairGeoMedia* media = geoFace->getMedia();
  static FairGeoBuilder* geoBuild = geoLoad->getGeoBuilder();
  
  FairGeoMedium* ShipMedium = media->getMedium(name);
  
  if (!ShipMedium) {
    Fatal("InitMedium", "Material %s not defined in media file.", name);
    return -1111;
  }
  TGeoMedium* medium = gGeoManager->GetMedium(name);
  if (medium != NULL)
    return ShipMedium->getMediumIndex();
  
  return geoBuild->createMedium(ShipMedium);
}

// -------------------------------------------------------------------------
/*
 *
 * @brief Processes a hit in the veto detector.
 *
 * This method is called whenever a hit is registered in the veto detector.
 * It processes the hit information and records the relevant data.
 *
 * @param x X-coordinate of the hit.
 * @param y Y-coordinate of the hit.
 * @param z Z-coordinate of the hit.
 * @param time Time of the hit.
 * @param energy Energy deposited by the hit.
 * @return True if the hit was processed successfully, false otherwise.
 */
Bool_t veto::ProcessHits(FairVolume* vol) {
  /*
   * This method is called from the MC stepping 
   */
  // Set parameters at entrance of volume. Reset ELoss.
  if (gMC->IsTrackEntering()) {
    fELoss = 0.;
    fTime = gMC->TrackTime() * 1.0e09;
    fLength = gMC->TrackLength();
    gMC->TrackPosition(fPos);
    gMC->TrackMomentum(fMom);
  }
  // Sum energy loss for all steps in the active volume
  fELoss += gMC->Edep();

  // Create vetoPoint at exit of active volume
  if (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()) {
    if (fELoss == 0.) {
      return kFALSE;
    }

    fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
    
    Int_t veto_uniqueId;
    gMC->CurrentVolID(veto_uniqueId);
    TParticle* p = gMC->GetStack()->GetCurrentTrack();
    Int_t pdgCode = p->GetPdgCode();
    TLorentzVector pos;
    gMC->TrackPosition(pos);
    TLorentzVector Mom;
    gMC->TrackMomentum(Mom);
    Double_t xmean = (fPos.X() + pos.X()) / 2.;
    Double_t ymean = (fPos.Y() + pos.Y()) / 2.;
    Double_t zmean = (fPos.Z() + pos.Z()) / 2.;
    AddHit(fTrackID,
	   veto_uniqueId,
	   TVector3(xmean, ymean, zmean),
	   TVector3(fMom.Px(), fMom.Py(), fMom.Pz()),
	   fTime,
	   fLength,
	   fELoss,
	   pdgCode,
	   TVector3(pos.X(), pos.Y(), pos.Z()),
	   TVector3(Mom.Px(), Mom.Py(), Mom.Pz()));
    
    // Increment number of veto det points in TParticle
    ShipStack* stack = (ShipStack*)gMC->GetStack();
    stack->AddPoint(kVETO);
  }

  return kTRUE;
}

void veto::EndOfEvent() {
  fvetoPointCollection->Clear();
}

void veto::PreTrack() {
  if (!fFastMuon) {
    return;
  }
  if (TMath::Abs(gMC->TrackPid()) != 13) {
    gMC->StopTrack();
  }
}

void veto::Register() {   // create a branch in the output tree
  FairRootManager::Instance()->Register("vetoPoint",
					"veto",
					fvetoPointCollection,
					kTRUE);   // kFALSE -> this collection will not be written to the file, will exist only during simulation.
}

TClonesArray* veto::GetCollection(Int_t iColl) const {
  if (iColl == 0) {
    return fvetoPointCollection;
  } else {
    return NULL;
  }
}

void veto::Reset() {
  fvetoPointCollection->Clear();
}

/*
 *
 * @brief Constructs the detector geometry.
 *
 * This function is responsible for setting up the geometry of the DecayVolume+SBT detector.
 * It is called during the detector's construction phase.
 */

void veto::ConstructGeometry() {
  TGeoVolume* top = gGeoManager->GetTopVolume();
  
  InitMedium("vacuum");
  InitMedium("helium");
  
  gGeoManager->SetNsegments(100);
  vetoMed = gGeoManager->GetMedium(vetoMed_name);   //! medium of veto counter, liquid or plastic scintillator
  decayVolumeMed = gGeoManager->GetMedium(decayVolumeMed_name);   // decay volume, air/helium/vacuum

  LOG(INFO) << "veto: Veto medium set as: " <<  vetoMed_name;
  LOG(INFO) << "veto: Decay Volume medium set as: " <<  decayVolumeMed_name;
  TGeoVolume* tDecayVol = new TGeoVolumeAssembly("DecayVolume");

  TGeoVolume* seg = MakeSegments();
  tDecayVol->AddNode(seg, 1, new TGeoTranslation(0, 0, 0));
  top->AddNode(tDecayVol, 1, new TGeoTranslation(0, 0, zStartDecayVol));   //));

  // only for fastMuon simulation, otherwise output becomes too big
  if (fFastMuon && fFollowMuon) {
    const char* Vol = "TGeoVolume";
    const char* Cavern = "Cavern";
    const char* Ain = "AbsorberAdd";
    const char* Aout = "AbsorberAddCore";
    TObjArray* volumelist = gGeoManager->GetListOfVolumes();
    int lastvolume = volumelist->GetLast();
    int volumeiterator = 0;
    while (volumeiterator <= lastvolume) {
      const char* volumename = volumelist->At(volumeiterator)->GetName();
      const char* classname = volumelist->At(volumeiterator)->ClassName();
      if (strstr(classname, Vol)) {
	if (strstr(volumename, Cavern) || strstr(volumename, Ain) || strstr(volumename, Aout)) {
	  AddSensitiveVolume(gGeoManager->GetVolume(volumename));
	  LOG(INFO) << "veto: made sensitive for following muons: " << volumename;
	}
      }
      volumeiterator++;
    }
  }
}

vetoPoint* veto::AddHit(Int_t trackID,
                        Int_t detID,
                        TVector3 pos,
                        TVector3 mom,
                        Double_t time,
                        Double_t length,
                        Double_t eLoss,
                        Int_t pdgCode,
                        TVector3 Lpos,
                        TVector3 Lmom) {
  TClonesArray& clref = *fvetoPointCollection;
  Int_t size = clref.GetEntriesFast();
  return new (clref[size]) vetoPoint(trackID, detID, pos, mom, time, length, eLoss, pdgCode, Lpos, Lmom);
}
