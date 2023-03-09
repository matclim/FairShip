#include "splitcalHit.h"
#include "splitcal.h"
#include "TVector3.h"
#include "FairRun.h"
#include "FairRunSim.h"
#include "TMath.h"
#include "TRandom1.h"
#include "TRandom3.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TGeoNavigator.h"
#include "TGeoShape.h"
#include "TGeoBBox.h"
#include "TF1.h"

#include <iostream>
#include <math.h>
#include <numeric>
#include <array>
using std::cout;
using std::endl;

Double_t speedOfLight = TMath::C() *100./1000000000.0 ; // from m/sec to cm/ns
// -----   Default constructor   -------------------------------------------
splitcalHit::splitcalHit()
  : ShipHit()
{
 flag = true;
}
// -----   Standard constructor   ------------------------------------------
splitcalHit::splitcalHit(Int_t detID, Float_t tdc)
  : ShipHit(detID,tdc)
{
 flag = true;
}

// -----   constructor from reconstructed positions, used for HPL extra hits   ------------------------------------------

splitcalHit::splitcalHit(Double_t recoX, Double_t recoY, Double_t recoZ, Double_t recoenergy, Double_t pex, Double_t pey)
  : ShipHit()
{
SetXYZ(recoX,recoY,recoZ);
	SetEnergy(recoenergy);
	SetisPrecisionLayer(1);
	px=pex;
	py=pey;
}

// -----   constructor from splitcalPoint   ------------------------------------------
splitcalHit::splitcalHit(splitcalPoint* p, Double_t t0)
  : ShipHit()
{

  flag = true;

  double pointX =  p->GetX();
  double pointY =  p->GetY();
  double pointZ =  p->GetZ();
  double pointT =  p->GetTime();
  double pointE =  p->GetEnergyLoss();
  TVector3 *pointPVector = new TVector3(p->GetPx(),p->GetPy(),p->GetPz());
  double pointPy = p->GetPy();
  double pointPz = p->GetPz();
  _detID =  p->GetDetectorID();

  //fdigi = t0 + t; 
  fdigi = t0 ;  
  // SetDigi(SetTimeRes(fdigi));
  SetDetectorID(_detID);

  TGeoNavigator* navigator = gGeoManager->GetCurrentNavigator();
  navigator->cd("cave/SplitCalDetector_1");
  TGeoVolume* caloVolume = navigator->GetCurrentVolume();
  // caloVolume->PrintNodes();

  int isPrec, nL, nMx, nMy, nS;
  std::string stripName = GetDetectorElementName(isPrec,_detID); // it also sets if strip gives x or y coordinate
  if(_detID>10000000){isPrec=1;}
  else isPrec=0;

  Decoder(_detID, nL, nMx, nMy, nS);

  
  SetIDs(isPrec, nL, nMx, nMy, nS);
  if(!isPrec){
  TGeoNode* strip = caloVolume->GetNode(stripName.c_str()); 

  const Double_t* stripCoordinatesLocal = strip->GetMatrix()->GetTranslation();
  Double_t stripCoordinatesMaster[3] = {0.,0.,0.};
  navigator->LocalToMaster(stripCoordinatesLocal, stripCoordinatesMaster);

  TGeoBBox* box = (TGeoBBox*)strip->GetVolume()->GetShape();
  double xHalfLength = box->GetDX();
  double yHalfLength = box->GetDY();
  double zHalfLength = box->GetDZ();
  
  TGeoNode* passiveLayer = caloVolume->GetNode("ECALfilter_200000"); // they are all the same
  TGeoBBox* boxPassive = (TGeoBBox*)passiveLayer->GetVolume()->GetShape();
  double zPassiveHalfLength = box->GetDZ();

  // std::cout<< "----------------------"<<std::endl;
  // std::cout<< "-- pointX = " << pointX << std::endl; 
  // std::cout<< "-- pointY = " << pointY << std::endl; 
  // std::cout<< "-- pointZ = " << pointZ << std::endl; 
  // std::cout<< "-- detID = " << detID << std::endl;
  // std::cout<< "-- stripName = " << stripName << std::endl;
  // std::cout<< "-- isPrec = " << isPrec << std::endl;
  // std::cout<< "-- nL = " << nL << std::endl;
  // std::cout<< "-- nMx = " << nMx << std::endl;
  // std::cout<< "-- nMy = " << nMy << std::endl;
  // std::cout<< "-- nS = " << nS << std::endl;    
  // std::cout<< "-- stripCoordinatesLocal[0] = " << stripCoordinatesLocal[0] << std::endl;
  // std::cout<< "-- stripCoordinatesLocal[1] = " << stripCoordinatesLocal[1] << std::endl;
  // std::cout<< "-- stripCoordinatesLocal[2] = " << stripCoordinatesLocal[2] << std::endl;
  // std::cout<< "-- stripCoordinatesMaster[0] = " << stripCoordinatesMaster[0] << std::endl;
  // std::cout<< "-- stripCoordinatesMaster[1] = " << stripCoordinatesMaster[1] << std::endl;
  // std::cout<< "-- stripCoordinatesMaster[2] = " << stripCoordinatesMaster[2] << std::endl;


  // TGeoNode* check = navigator->FindNode(pointX,pointY,pointZ);

  SetEnergy(pointE);
  SetXYZ(stripCoordinatesMaster[0], stripCoordinatesMaster[1], stripCoordinatesMaster[2]);
  SetXYZErrors(xHalfLength,yHalfLength,2*(zHalfLength+zPassiveHalfLength));
 }
 else{
  ComputeCloudSigma(pointE);
  ComputeHPL(pointX,pointY,pointZ,pointE,pointPVector,_detID);

	} 

 delete pointPVector;
}




void splitcalHit::ComputeCloudSigma(double pointE){

// If there is energy dependence, it should be implemented
// Values based off of 2019 prototype MicroMegas
// 10% threshold 0.00170933 GeV (purely arbitrary)
//
//
// drift region: ~1.1cm, E~600 V/cm (not used, cloud assumed insignificant)
// electron speed 1.4 * 10^7 m/s
// drift time 360 ps
// sigma = 5.4 microm
//
// amplification region ~100 microm + 50 for Y, E~50 kV/cm
// electron speed 4*10^7 m/s (42163702 m/s)
// drift time 25 ps (2.5*10^-11 s) / 36 ps for Y
// D Ar = 0.04 cm^2/s (no interactions assumed with 7% of CO2)
// sigma = sqrt(2Dt) = 1.4 microm
// sigmaY = 1.7 microm
// For now just sum the cloud widths, a better system is needed though

	_cloudsigmaX = 0.00014 + 0.00054;
	_cloudsigmaY = 0.00017 + 0.00054;
	
}


Double_t splitcalHit::ComputeHPLDiffusedEnergy(double x1, double x2, double y1, double y2, double mux, double muy){
//x1 & y1: lower limit 
//x2 & y2 : upper limit
//mux & muy real center of the gaussian	

	Double_t e_fraction(0.25*TMath::Abs(TMath::Erf((x1-mux)/_cloudsigmaX) - TMath::Erf((x2-mux)/_cloudsigmaX))*TMath::Abs(TMath::Erf((y1-muy)/_cloudsigmaY) - TMath::Erf((y2-muy)/_cloudsigmaY)));

	return e_fraction;


}


//average energy in HPL 0.00086
//Convert linearily to 10 bit, arbitrary fix of average at 400
//a = 465116
std::vector<double> splitcalHit::DigitizeEnergy(std::vector<double> vecofe){
	std::vector<double> vecofdig;
	for(int i=0;i<vecofe.size();i++){
		vecofdig.push_back(vecofe[i]*465116);
	}
	return vecofdig;	

}



std::array<double,2> splitcalHit::ComputeBaryCenter(std::vector<std::array<double,4>> vecofrec,std::vector<double> vecofe){
	double sumX(0);
	double sumY(0);
	double EXsum(0);
	double EYsum(0);
	for(int i=0;i<vecofrec.size();++i){
		sumX+=vecofrec[0][i]*vecofe[i];
		sumY+=vecofrec[1][i]*vecofe[i];
		EXsum+=vecofe[i];
		EYsum+=vecofe[i];
	}
	std::array<double,2> recoposi={sumX/EXsum,sumY/EYsum};
	return recoposi;
}





//#same as now generate ionisations
//rng for number of ionisations on the track
//rng for energy proportion
//then project on x&y with amplification
//then (later) add clouds and diffusion

void splitcalHit::SelectHPLayer(int DETID){

	switch(DETID){

	case 100700000: HPLayer =0; break;
	case 101100000: HPLayer =1; break;
	case 101400000: HPLayer =2; break;
	default: HPLayer=-1;
	}


}



void splitcalHit::ComputeHPL(double pointX, double pointY, double pointZ, double pointE, TVector3 *pointPVector, int DETID){
	px=pointX;
	py=pointY;

	SelectHPLayer(DETID);

        double borderenergy[4]; //energy in bordering strips because of diffusion 0= left strip 1= right strip 2= lower strip 3 = top strip
	double recX,recY,recZ; //reconstructed positions
	std::vector<std::array<double,4>> vecofrec; // vector of strip reconstruction
	std::vector<std::array<double,3>> vecofspot; //vector of different strip positions
	TRandom3 *rgen = new TRandom3();	
	std::vector<std::array<double,3>> vecofpos; //vector of point derived position after drift
	std::vector<double> vecofrand;	           // vector of random energy proportion and energy
	int ion_no = rgen->Poisson(int(round(33))); //ionization number average number of primary interactions is 29.4/cm
	double ion_point=0;			    // position on the drift where the interaction takes place
	for(int ioniz=0;ioniz<ion_no;ioniz++){
		//position simulation

		ion_point=rgen->Rndm()*pointPVector->Mag(); //
		TVector3 *shortvector= new TVector3(*pointPVector);
		shortvector->SetMag(ion_point);
		vecofpos.push_back({pointX + shortvector->Px()*1.12,pointY + shortvector->Py()*1.12,pointZ});
		
		//energy simulation step 1
		vecofrand.push_back(rgen->Rndm());
		delete shortvector;		
		}
		
		float totalrand=0;
		for(auto ind : vecofrand){
			totalrand+=ind;
		}

		for(auto &ind : vecofrand){
			ind=ind/totalrand;
		}
		for(int i=0;i<ion_no;i++){
			if(TMath::Abs(vecofpos[i][0])<600 && TMath::Abs(vecofpos[i][1]) <1200){
				recX=round(vecofpos[i][0]/0.05);		
				recY=round(vecofpos[i][1]/0.05);		
				recX= recX * 0.05;
				recY= recY * 0.05;
				vecofrec.push_back({recX,recY,pointE*vecofrand[i]*ComputeHPLDiffusedEnergy(recX-0.025,recX+0.025,recY-0.025,recY+0.025,vecofpos[i][0],vecofpos[i][1]),0});

				borderenergy[0]=pointE*vecofrand[i]*ComputeHPLDiffusedEnergy(recX-0.075,recX-0.025,recY-0.025,recY+0.025,vecofpos[i][0],vecofpos[i][1]);
				borderenergy[1]=pointE*vecofrand[i]*ComputeHPLDiffusedEnergy(recX+0.025,recX+0.075,recY-0.025,recY+0.025,vecofpos[i][0],vecofpos[i][1]);
				borderenergy[2]=pointE*vecofrand[i]*ComputeHPLDiffusedEnergy(recX-0.025,recX+0.025,recY-0.075,recY-0.025,vecofpos[i][0],vecofpos[i][1]);
				borderenergy[3]=pointE*vecofrand[i]*ComputeHPLDiffusedEnergy(recX-0.025,recX+0.025,recY+0.025,recY+0.075,vecofpos[i][0],vecofpos[i][1]);
				if(borderenergy[0]>HPLentresh){
					vecofrec.push_back({recX-0.05,recY,borderenergy[0],0});
					sizeofExtraDiffusionHits+=1;
				}			
				if(borderenergy[1]>HPLentresh){
					vecofrec.push_back({recX+0.05,recY,borderenergy[1],0});
					sizeofExtraDiffusionHits+=1;
				}			
				if(borderenergy[2]>HPLentresh){
					vecofrec.push_back({recX,recY-0.05,borderenergy[2],0});
					sizeofExtraDiffusionHits+=1;
				}			
				if(borderenergy[3]>HPLentresh){
					vecofrec.push_back({recX,recY+0.05,borderenergy[3],0});
					sizeofExtraDiffusionHits +=1;
				}			


//				vecofrec.push_back({recX,recY,pointE*vecofrand[i],0});
				}
		}

		//selection loop
		std::vector<int> spotters;
		for(int i=0;i<vecofrec.size();i++){ //the i loop sets the reference, 0,1 are x,y ; 2 is energy ; 3 is the status [0=unused, 1=used, 2= used and ref]
			if(vecofrec[i][3]==0){ 
				vecofrec[i][3]=2;
				spotters.push_back(i);
				for(int j=0;j<vecofrec.size();j++){
					if(vecofrec[j][3]==0 && vecofrec[i][0]==vecofrec[j][0] && vecofrec[i][1]==vecofrec[j][1]){
						vecofrec[i][2] += vecofrec[j][2];	
						vecofrec[j][3]=1;
					}
				}
			}	
		}

		//with barycenter reconstruction
		std::vector<double> vecoffinalenergy;
		std::vector<std::array<double,2>> vecoftpos;

		for(int i=0;i<vecofrec.size();i++){
			vecoffinalenergy.push_back(vecofrec[i][2]);
			vecoftpos.push_back({vecofrec[i][0],vecofrec[i][1]});
		}

			

		//without barycenter reconstruction
		SetXYZ(vecofrec[0][0],vecofrec[0][1],pointZ);
		SetEnergy(vecofrec[0][2]);
		firstextraflag=kTRUE;
		//cout << "Energy " << vecofrec[0][2] << " pointE " << pointE << endl;
		for(int i=1;i<vecofrec.size();i++){	
			if(vecofrec[i][3]==2){
				//cout << "i = " << i << " " << endl;
				//splitcalHit *b= new splitcalHit(vecofrec[i][0],vecofrec[i][1],pointZ,vecofrec[i][2],);
				_vecofExtraHits.push_back(new splitcalHit(vecofrec[i][0],vecofrec[i][1],pointZ,vecofrec[i][2],pointX,pointY));
				//delete b;
			}
		}
		sizeofExtraHits	= _vecofExtraHits.size();

	delete rgen;
}

void splitcalHit::ComputeHPLInterim(double pointX, double pointY, double pointZ){
	double recX,recY,recZ;







	if(TMath::Abs(pointX)<600 && TMath::Abs(pointY) <1200){
		recX = round(pointX/0.00029);
		if(recX-pointX)
		recY = round(pointY/0.000475);
		recX = recX * 0.00029;
		recY = recY * 0.000475;}
	else{recX=600; recY=1200;}

	SetXYZ(recX,recY,pointZ);
	SetXYZErrors(TMath::Abs(recX-pointX),TMath::Abs(recY-pointY),0);

}

std::string splitcalHit::GetPaddedString(int& id){ 

  //zero padded string 
  int totalLength = 9; 
  std::string stringID = std::to_string(id);
  std::string encodedID = std::string(totalLength - stringID.length(), '0') + stringID;
 
 return encodedID;

}

std::string splitcalHit::GetDetectorElementName(int isPrec,int& id){

  std::string encodedID = GetPaddedString(id);
  //std::cout << "-- encodedID = " << encodedID <<std::endl;
  int nL, nMx, nMy, nS;
  Decoder(encodedID, nL, nMx, nMy, nS);

  std::string name;
  if (isPrec==1) {
    name = "ECALdet_gas_";
    SetIsX(true);
    SetIsY(true);
  } else if (nL%2==0) {
    name = "stripGivingY_"; 
    SetIsX(false);
    SetIsY(true);
  } else {
    name = "stripGivingX_";
    SetIsX(true);
    SetIsY(false);
  }
  name = name + std::to_string(id);
  // std::cout << "--GetDetectorElementName - name  = " << name <<std::endl;

  return name;

}

void splitcalHit::Decoder(std::string& encodedID, int& nLayer, int& nModuleX,  int& nModuleY, int& nStrip){

  std::string subtring;

  //subtring = encodedID.substr(0, 1);
  //isPrecision = atoi(subtring.c_str());

  subtring = encodedID.substr(1,3);
  nLayer = atoi(subtring.c_str());

  subtring = encodedID.substr(4,1); 
  nModuleX = atoi(subtring.c_str());

  subtring = encodedID.substr(5,1); 
  nModuleY = atoi(subtring.c_str());

  subtring = encodedID.substr(6,3);
  nStrip = atoi(subtring.c_str());

}

void splitcalHit::Decoder(int& id, int& nLayer, int& nModuleX,  int& nModuleY, int& nStrip){

  std::string encodedID = GetPaddedString(id);
  Decoder(encodedID, nLayer, nModuleX, nModuleY, nStrip);

}


double splitcalHit::GetEnergyWeightForIndex(int index){

  int iw = 0;
  for(size_t i=0; i<_vecClusterIndices.size(); i++) {
    if (_vecClusterIndices.at(i) == index) {
      iw = i;
      break; 
    }
  }
  return _vecEnergyWeights.at(iw); 
}


double splitcalHit::GetEnergyForCluster(int i){

  double unweightedEnergy = GetEnergy();
  double weight = GetEnergyWeightForIndex(i);
  double energy = unweightedEnergy*weight;
  return energy;

}

// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
splitcalHit::~splitcalHit() { 
	for(auto &a : _vecofExtraHits){
		delete a;
	}

}
// -------------------------------------------------------------------------

// -----   Public method Print   -------------------------------------------
void splitcalHit::Print() const
{
  //  cout << "-I- splitcalHit: splitcal hit " << " in detector " << fDetectorID << endl;
  //  cout << "  TDC " << fdigi << " ns" << endl;
  // std::cout<< "-I- splitcalHit: " <<std::endl;
  std::cout<< "------- " <<std::endl;
  std::cout<< "    (x,y,z) = " 
	   << _x << " +- " << _xError << " ,  "
	   << _y << " +- " << _yError << " ,  "
	   << _z << " +- " << _zError <<std::endl;  
  std::cout<< "    isP, nL, nMx, nMy, nS = " 
	   << _isPrecisionLayer << " , "  
	   << _nLayer << " , "  
	   << _nModuleX << " , " 
	   << _nModuleY << " , " 
	   << _nStrip << std::endl;
   std::cout<< "------- " <<std::endl;
}
// -------------------------------------------------------------------------

ClassImp(splitcalHit)

