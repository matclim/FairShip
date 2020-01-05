#ifndef PIXELMODULES_H
#define PIXELMODULES_H

#include "FairModule.h"                 // for FairModule
#include "FairDetector.h"                  // for FairDetector

#include "Rtypes.h"                     // for ShipMuonShield::Class, Bool_t, etc

#include <string>                       // for string

#include "TVector3.h"
#include "TLorentzVector.h"
#include "unordered_map"

class PixelModulesPoint;
class FairVolume;
class TClonesArray;

class PixelModules:public FairDetector
{
  public:
  PixelModules(const char* name, const Double_t DX, const Double_t DY, const Double_t DZ,Bool_t Active, const char* Title="PixelModules", const Int_t nSl=1);
    PixelModules();
    virtual ~PixelModules();
    void ConstructGeometry();
    void SetZsize(const Double_t MSsize);
    void SetBoxParam(Double_t SX, Double_t SY, Double_t SZ, Double_t zBox, Double_t SZPixel, Double_t D1short, Double_t D1long,Double_t SiliconDZthin, Double_t SiliconDZthick,Int_t nstation,Double_t posx,Double_t posy,Double_t posz,Double_t anglex,Double_t angley,Double_t anglez,Int_t nSl=1);
//    SetBoxParam(DX,DY,DZ, zBox, DimZpixelbox, D1short, D1long,DimZSithin, DimZSithick,nSlice)
void SetBoxParam(Double_t SX, Double_t SY, Double_t SZ, Double_t zBox, Double_t SZPixel, Double_t D1short, Double_t D1long,Double_t SiliconDZthin, Double_t SiliconDZthick,Int_t nstation,Double_t posx,Double_t posy,Double_t posz,Double_t anglex,Double_t angley,Double_t anglez);
    
    void SetBoxParam(Double_t SX, Double_t SY, Double_t SZ, Double_t zBox,Double_t SZPixel, Double_t Dim1Short, Double_t Dim1Long);
    void SetSiliconDZ(Double_t SiliconDZthin, Double_t SiliconDZthick);  
    void SetSiliconStationPositions(Int_t nstation, Double_t posx, Double_t posy, Double_t posz);
    void SetSiliconStationAngles(Int_t nstation, Double_t anglex, Double_t angley, Double_t anglez);
    void SetSiliconDetNumber(Int_t nSilicon);
    void SetSiliconSlicesNumber(Int_t nSl=1);
    void ComputeDimZSlice();
    /**      Initialization of the detector is done here    */
    virtual void Initialize();

    /**       this method is called for each step during simulation
     *       (see FairMCApplication::Stepping())
     */
    virtual Bool_t ProcessHits( FairVolume* v=0);

    /**       Registers the produced collections in FAIRRootManager.     */
    virtual void   Register();

    /** Gets the produced collections */
    virtual TClonesArray* GetCollection(Int_t iColl) const ;

    /**      has to be called after each event to reset the containers      */
    virtual void   Reset();

    /**      This method is an example of how to add your own point
     *       of type muonPoint to the clones array
     */
    PixelModulesPoint* AddHit(Int_t trackID, Int_t detID,
                         TVector3 pos, TVector3 mom,
                         Double_t time, Double_t length,
                         Double_t eLoss, Int_t pdgCode);

    /** The following methods can be implemented if you need to make
     *  any optional action in your detector during the transport.
     */

    virtual void   CopyClones( TClonesArray* cl1,  TClonesArray* cl2 ,
                              Int_t offset) {;}
    virtual void   SetSpecialPhysicsCuts() {;}
    virtual void   EndOfEvent();
    virtual void   FinishPrimary() {;}
    virtual void   FinishRun() {;}
    virtual void   BeginPrimary() {;}
    virtual void   PostTrack() {;}
    virtual void   PreTrack() {;}
    virtual void   BeginEvent() {;}

    void DecodeVolumeID(Int_t detID,int &nHPT);

private:

    /** Track information to be stored until the track leaves the
     active volume.
     */
    Int_t          fTrackID;           //!  track index
    Int_t          fPdgCode;           //!  pdg code
    Int_t          fVolumeID;          //!  volume id
    TLorentzVector fPos;               //!  position at entrance
    TLorentzVector fMom;               //!  momentum at entrance
    Double32_t     fTime;              //!  time
    Double32_t     fLength;            //!  length
    Double32_t     fELoss;             //!  energy loss
    /*Number of Silicon pixel segments, 10 is the recommended amount for digitized runs, 1 for simple runs, to be set in charm geometry file*/
    Int_t nSlices=10; 
    
    /** container for data points */
    TClonesArray*  fPixelModulesPointCollection;

    Int_t InitMedium(const char* name);
    Int_t *GetIDlist1(); 
    Int_t *GetIDlist10();
    Bool_t *vertical; 

protected:

    Double_t *GetPosSize1();
    Double_t *GetPosSize10();
    void SetPositionSize();
    Double_t Dim1Short, Dim1Long;
    Int_t *SetIDs();   
    void SetVertical(); 
    Bool_t *GetVertical1();
    Bool_t *GetVertical10();
 
    Double_t SBoxX = 0;
    Double_t SBoxY = 0;
    Double_t SBoxZ = 0;

    Double_t BoxX = 0;
    Double_t BoxY = 0;
    Double_t BoxZ = 0;
    Double_t zBoxPosition = 0;

    Double_t DimX =0;
    Double_t DimY =0;
    Double_t DimZ = 0;
    Double_t DimZWindow=0.0110;
    Double_t Windowx=5;
    Double_t Windowy=5;
    Double_t FrontEndthick=0.0150;
    Double_t FlexCuthick=0.0100;
    Double_t FlexKapthick=0.0050;
    Double_t overlap=0;
    Double_t DimZPixelBox;
    Int_t *PixelIDlist;

    static const Int_t nSi1=12;
    static const Int_t nSi10=120;
    Int_t numSi;
    Double_t DimZSithin=0.02;
    Double_t DimZSithick=0.0245;
    Double_t DimZThinSlice;
    Double_t DimZThickSlice;
    Double_t *xs, *ys, *zs;
    Double_t *xangle, *yangle, *zangle;
    
    PixelModules(const PixelModules&);
    PixelModules& operator=(const PixelModules&);
    ClassDef(PixelModules,1)

};
#endif
