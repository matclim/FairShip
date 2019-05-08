// PixelModules.cxx
//  PixelModules, twelve pixel modules physically connected two by two.

#include "PixelModules.h"
//#include "MagneticPixelModules.h" 
#include "PixelModulesPoint.h"
#include "TGeoManager.h"
#include "FairRun.h"                    // for FairRun
#include "FairRuntimeDb.h"              // for FairRuntimeDb
#include <iosfwd>                    // for ostream
#include "TList.h"                      // for TListIter, TList (ptr only)
#include "TObjArray.h"                  // for TObjArray
#include "TString.h"                    // for TString
#include "TClonesArray.h"
#include "TVirtualMC.h"

#include "TGeoBBox.h"
#include "TGeoTrd1.h"
#include "TGeoCompositeShape.h"
#include "TGeoTube.h"
#include "TGeoArb8.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TParticle.h"
#include "TVector3.h"

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
#include "ShipUnit.h"
#include "ShipStack.h"

#include "TGeoUniformMagField.h"
#include <stddef.h>                     // for NULL
#include <iostream>                     // for operator<<, basic_ostream, etc

using std::cout;
using std::endl;
using namespace ShipUnit;

PixelModules::PixelModules()
  : FairDetector("HighPrecisionTrackers",kTRUE, kPixelModules),
    numSi(nSi),
    fTrackID(-1),
    fPdgCode(),
    fVolumeID(-1),
    fPos(),
    fMom(),
    fTime(-1.),
    fLength(-1.),
    fELoss(-1),
    fPixelModulesPointCollection(new TClonesArray("PixelModulesPoint"))
{
}

PixelModules::PixelModules(const char* name, const Double_t DX, const Double_t DY, const Double_t DZ, Bool_t Active,const char* Title)
  : FairDetector(name, Active, kPixelModules),
    numSi(nSi),
    fTrackID(-1),
    fPdgCode(),
    fVolumeID(-1),
    fPos(),
    fMom(),
    fTime(-1.),
    fLength(-1.),
    fELoss(-1),
    fPixelModulesPointCollection(new TClonesArray("PixelModulesPoint"))
{ 
  DimX = DX;
  DimY = DY;
  DimZ = DZ;
}

PixelModules::~PixelModules()
{
    if (fPixelModulesPointCollection) {
        fPixelModulesPointCollection->Delete();
        delete fPixelModulesPointCollection;
    }
}

void PixelModules::Initialize()
{
    FairDetector::Initialize();
}

// -----   Private method InitMedium 
Int_t PixelModules::InitMedium(const char* name)
{
   static FairGeoLoader *geoLoad=FairGeoLoader::Instance();
   static FairGeoInterface *geoFace=geoLoad->getGeoInterface();
   static FairGeoMedia *media=geoFace->getMedia();
   static FairGeoBuilder *geoBuild=geoLoad->getGeoBuilder();

   FairGeoMedium *ShipMedium=media->getMedium(name);

   if (!ShipMedium)
   {
     Fatal("InitMedium","Material %s not defined in media file.", name);
     return -1111;
   }
   TGeoMedium* medium=gGeoManager->GetMedium(name);
   if (medium!=NULL)
     return ShipMedium->getMediumIndex();
   return geoBuild->createMedium(ShipMedium);
}

void PixelModules::SetBoxParam(Double_t SX, Double_t SY, Double_t SZ, Double_t zBox, Double_t SZPixel, Double_t D1short, Double_t D1long)
{
  SBoxX = SX;
  SBoxY = SY;
  SBoxZ = SZ;
  zBoxPosition = zBox;
  DimZPixelBox = SZPixel;
  Dim1Short = D1short;
  Dim1Long = D1long;
}

void PixelModules::SetSiliconDZ(Double_t SiliconDZthin,Double_t SiliconDZthick)
{
  DimZSithin = SiliconDZthin;
  DimZSithick= SliconDZthick
}


void PixelModules::SetSiliconStationPositions(Int_t nstation, Double_t posx, Double_t posy, Double_t posz)
{
 xs[nstation] = posx;
 ys[nstation] = posy;
 zs[nstation] = posz;
}

void PixelModules::SetSiliconStationAngles(Int_t nstation, Double_t anglex, Double_t angley, Double_t anglez)
{
 xangle[nstation] = anglex;
 yangle[nstation] = angley;
 zangle[nstation] = anglez;
}

void PixelModules::SetSiliconSlicesNumber(Int_t nSl)
{
nSlices=nSl;
}

void PixelModules::ComputeDimZSlice(){
DimZSlice=DimZSi/nSlices;
}


//
void PixelModules::ConstructGeometry()
{ 
    InitMedium("air");
  TGeoMedium *air = gGeoManager->GetMedium("air");

    InitMedium("iron");
    TGeoMedium *Fe =gGeoManager->GetMedium("iron");
    
    InitMedium("silicon");
    TGeoMedium *Silicon = gGeoManager->GetMedium("silicon");

    InitMedium("aluminium");
    TGeoMedium *Aluminium = gGeoManager->GetMedium("aluminium");

    InitMedium("CoilCopper");
    TGeoMedium *Cu  = gGeoManager->GetMedium("CoilCopper");

    InitMedium("CoilAluminium");
    TGeoMedium *Al  = gGeoManager->GetMedium("CoilAluminium");

    InitMedium("TTmedium");
    TGeoMedium *TT  = gGeoManager->GetMedium("TTmedium");
    
    InitMedium("STTmix8020_2bar");
    TGeoMedium *sttmix8020_2bar   = gGeoManager->GetMedium("STTmix8020_2bar");
  
    TGeoVolume *top = gGeoManager->GetTopVolume();

    //computing the largest offsets in order to set PixelBox dimensions correctly
    Double_t offsetxmax = 0., offsetymax = 0.;
    for (int istation = 0; istation < nSi; istation++){
     if (TMath::Abs(xs[istation]) > offsetxmax) offsetxmax = TMath::Abs(xs[istation]);
     if (TMath::Abs(ys[istation]) > offsetymax) offsetymax = TMath::Abs(ys[istation]);
    }
    TGeoBBox *PixelBox = new TGeoBBox("PixelBox", Dim1Long/2 + offsetxmax, Dim1Long/2 + offsetymax, DimZPixelBox/2.); //The box is symmetric, offsets are not. So we enlarge the offset by a factor two for coverage
    TGeoVolume *volPixelBox = new TGeoVolume("volPixelBox",PixelBox,air);
    Double_t inimodZoffset(zs[0]) ;//initial Z offset of Pixel Module 0 so as to avoid volume extrusion
    top->AddNode(volPixelBox, 1, new TGeoTranslation(0,0,zBoxPosition+ inimodZoffset)); //volume moved in
    

    TGeoBBox *Pixelythin = new TGeoBBox("Pixelythin", Dim1Short/2, Dim1Long/2, DimZThinSlice/2); //long along y
    TGeoVolume *volPixelythin = new TGeoVolume("volPixelythin",Pixelythin,Silicon); 
    volPixely->SetLineColor(kBlue-5);
    AddSensitiveVolume(volPixelythin);

    TGeoBBox *Pixelxthin = new TGeoBBox("Pixelx", (Dim1Long)/2, (Dim1Short)/2, DimZThinSlice/2); //long along x
    TGeoVolume *volPixelxthin = new TGeoVolume("volPixelxthin",Pixelxthin,Silicon); 
    volPixelx->SetLineColor(kBlue-5);
    AddSensitiveVolume(volPixelxthin);

    TGeoBBox *Pixelythick = new TGeoBBox("Pixelythick", Dim1Short/2, Dim1Long/2, DimZThinSlice/2); //long along y
    TGeoVolume *volPixelythick = new TGeoVolume("volPixelythick",Pixelythick,Silicon); 
    volPixely->SetLineColor(kBlue-5);
    AddSensitiveVolume(volPixelythick);

    TGeoBBox *Pixelxthick = new TGeoBBox("Pixelx", (Dim1Long)/2, (Dim1Short)/2, DimZThinSlice/2); //long along x
    TGeoVolume *volPixelxthick = new TGeoVolume("volPixelxthick",Pixelxthick,Silicon); 
    volPixelx->SetLineColor(kBlue-5);
    AddSensitiveVolume(volPixelxthick);
  
  ///////////////////////////////////////////////////////Passive material///////////////////////////////////////////////////////
   
    TGeoBBox *WindowBox = new TGeoBBox("WindowBox",Windowx/2, Windowy/2,DimZWindow/2);
    TGeoVolume *volWindow = new TGeoVolume("volWindow",WindowBox,Al);
    volWindow->SetLineColor(kGray);
    //AddSensitiveVolume(volWindow);

    TGeoBBox *PixelFramesx = new TGeoBBox("PixelFramesx",Dim1Short/2,Dim1Long/2, DimZAlu/2);
    TGeoVolume *volAluFramex = new TGeoVolume("volAluFramex",PixelFramesx,Aluminium);
    volAluFramex->SetLineColor(kGray);
   // AddSensitiveVolume(volAluFramex);

    TGeoBBox *PixelFramesy = new TGeoBBox("PixelFramesy",Dim1Short/2,Dim1Long/2, DimZAlu/2);
    TGeoVolume *volAluFramey = new TGeoVolume("volAluFramey",PixelFramesy,Aluminium);
    volAluFramey->SetLineColor(kGray);
    //AddSensitiveVolume(volAluFramey);

////////////////////////////////////////////////////////End passive material////////////////////////////////////////////////////////////////

    //id convention: 1{a}{b}{c}, a = number of pair (from 1 to 6), b = element of the pair (1 or 2)
    int chi=0;
    Int_t PixelIDlist[nSi];
    for(int i=1110;i<1130;i++){
    PixelIDlist[chi]=i;
    chi++;
    }
    
    for(int i=1210;i<1230;i++){
    PixelIDlist[chi]=i;
    chi++;
    }

    for(int i=1310;i<1330;i++){
    PixelIDlist[chi]=i;
    chi++;
    }
    for(int i=1410;i<1430;i++){
    PixelIDlist[chi]=i;
    chi++;
    }
    for(int i=1510;i<1530;i++){
    PixelIDlist[chi]=i;
    chi++;
    }
    for(int i=1610;i<1630;i++){
    PixelIDlist[chi]=i;
    chi++;
    }
    //Alternated pixel stations optimized for y and x measurements
    Bool_t vertical[nSi];
    for(int i=0;i<nSi;i++){
	vertical[i]=kFALSE;
	}
    for(int i=0;i<3;i++){
    	for(int j =0;j<20;j++){
	vertical[i*40+j]=kTRUE;
	}
    }
  volWindow->AddNode(0,0,new TGeoTranslation(0,0,-DimZPixelBox/2.-inimodZoffset));

    for (int ipixel = 0; ipixel < nSi; ipixel++){
      if (vertical[ipixel]){
		if(PixelIDlist[ipixel]) {
			if(PixelIDlist[ipixel]<1130 || (PixelIDlist[ipixel]>1219 && PixelIDlist[ipixel]<1620)) volPixelBox->AddNode(volPixelythick, PixelIDlist[ipixel], new TGeoTranslation(xs[ipixel],ys[ipixel],-DimZPixelBox/2.+ zs[ipixel]-inimodZoffset)); //compensation for the Node offset
			else volPixelBox->AddNode(volPixelythin, PixelIDlist[ipixel], new TGeoTranslation(xs[ipixel],ys[ipixel],-DimZPixelBox/2.+ zs[ipixel]-inimodZoffset)); //compensation for the Node offset
		else volPixelBox->AddNode(volPixely, 9000, new TGeoTranslation(xs[ipixel],ys[ipixel],-DimZPixelBox/2.+ zs[ipixel]-inimodZoffset)); //else here used for debugging, if the number of slices isn't 10

		if((ipixel+nSlices)%9==1) volAluFramex->AddNode(volAluFramex, 0,new TGeoTranslation(xs[ipixel],ys[ipixel],-DimZPixelBox/2.+ zs[ipixel]-inimodZoffset+DimZSlice));
		}
      else{ 
		if(PixelIDlist[ipixel]) 
			if(PixelIDlist[ipixel]<1130 || (PixelIDlist[ipixel]>1219 && PixelIDlist[ipixel]<1620))volPixelBox->AddNode(volPixelxthick, PixelIDlist[ipixel], new TGeoTranslation(xs[ipixel],ys[ipixel],-DimZPixelBox/2.+ zs[ipixel]-inimodZoffset));
			else volPixelBox->AddNode(volPixelxthin, PixelIDlist[ipixel], new TGeoTranslation(xs[ipixel],ys[ipixel],-DimZPixelBox/2.+ zs[ipixel]-inimodZoffset));
		
		else volPixelBox->AddNode(volPixelx, 9000, new TGeoTranslation(xs[ipixel],ys[ipixel],-DimZPixelBox/2.+ zs[ipixel]-inimodZoffset)); //debugging else for if nSlice!=10


		if((ipixel+nSlices)%9==1)volAluFramey->AddNode(volAluFramey, 0,new TGeoTranslation(xs[ipixel],ys[ipixel],-DimZPixelBox/2.+ zs[ipixel]-inimodZoffset+DimZSlice));
		}
	
	}
}



Bool_t  PixelModules::ProcessHits(FairVolume* vol)
{
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
    
    // Create muonPoint at exit of active volume
    if ( gMC->IsTrackExiting()    ||
        gMC->IsTrackStop()       ||
        gMC->IsTrackDisappeared()   ) {
        fTrackID  = gMC->GetStack()->GetCurrentTrackNumber();

        if (fELoss == 0. ) { return kFALSE; }
        TParticle* p=gMC->GetStack()->GetCurrentTrack();
        Int_t pdgCode = p->GetPdgCode();
	//Int_t fMotherID =p->GetFirstMother();
        gMC->CurrentVolID(fVolumeID);	

        TLorentzVector Pos; 
        gMC->TrackPosition(Pos); 
        Double_t xmean = (fPos.X()+Pos.X())/2. ;      
        Double_t ymean = (fPos.Y()+Pos.Y())/2. ;      
        Double_t zmean = (fPos.Z()+Pos.Z())/2. ;     

	AddHit(fTrackID, fVolumeID, TVector3(xmean, ymean,  zmean), TVector3(fMom.Px(), fMom.Py(), fMom.Pz()), fTime, fLength,fELoss, pdgCode);
        
        // Increment number of muon det points in TParticle
        ShipStack* stack = (ShipStack*) gMC->GetStack();
        stack->AddPoint(kPixelModules);
    }
    
    return kTRUE;
}

void PixelModules::EndOfEvent()
{
    fPixelModulesPointCollection->Clear();
}


void PixelModules::Register()
{
    
    /** This will create a branch in the output tree called
     PixelModulesPoint, setting the last parameter to kFALSE means:
     this collection will not be written to the file, it will exist
     only during the simulation.
     */
    
    FairRootManager::Instance()->Register("PixelModulesPoint", "PixelModules",
                                          fPixelModulesPointCollection, kTRUE);
}

// -----   Public method to Decode volume info  -------------------------------------------
// -----   returns hpt, arm, rpc numbers -----------------------------------
void PixelModules::DecodeVolumeID(Int_t detID,int &nHPT)
{
  nHPT = detID;
}

TClonesArray* PixelModules::GetCollection(Int_t iColl) const
{
    if (iColl == 0) { return fPixelModulesPointCollection; }
    else { return NULL; }
}

void PixelModules::Reset()
{
    fPixelModulesPointCollection->Clear();
}


PixelModulesPoint* PixelModules::AddHit(Int_t trackID, Int_t detID,
                        TVector3 pos, TVector3 mom,
                        Double_t time, Double_t length,
					    Double_t eLoss, Int_t pdgCode)

{
    TClonesArray& clref = *fPixelModulesPointCollection;
    Int_t size = clref.GetEntriesFast();

    return new(clref[size]) PixelModulesPoint(trackID, detID, pos, mom,time, length, eLoss, pdgCode);
}


ClassImp(PixelModules)
