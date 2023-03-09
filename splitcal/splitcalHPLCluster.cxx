#include "splitcalHPLCluster.h"
#include "TMath.h"

#include <iostream>
#include <math.h>
#include <functional>   
#include <numeric>
#include <map>


// -----   constructor from list/vector of splitcalHit   ------------------------------------------
// splitcalHPLCluster::splitcalHPLCluster(boost::python::list& l)
// {
//   std::vector<splitcalHit > v;
//   for(int i=0; i<boost::python::len(l); i++) { 
//     v.push_back(boost::python::extract<splitcalHit >(l[i]));
//   }
//   SetVectorOfHits(v);
// }

// -----   Default constructor   -------------------------------------------
splitcalHPLCluster::splitcalHPLCluster()
{
}
// -----   constructor from splitcalHit   ------------------------------------------
splitcalHPLCluster::splitcalHPLCluster(splitcalHit* h)
{
  _vectorOfHits.push_back(h);
 
 
}



double splitcalHPLCluster::GetX(int i){

	return _vectorOfHits.at(i)->GetX();


}

int splitcalHPLCluster::DigitizeEnergy(double energy){

	//digitize over 6 bits, assume 20 corresponds to average value of 0.00086 GeV
	return 23256*energy;
}



void splitcalHPLCluster::Analyze(){
	Clustersize = vectorOfClusters.size();
	double Xtopsum=0;
	double Ytopsum=0;
	double Zbary = vectorOfClusters.at(0).first->GetZ();
	Layer = vectorOfClusters.at(0).first->GetLayer();
	int energysum=0;
	for(auto itr=vectorOfClusters.begin();itr<vectorOfClusters.end();itr++){
		itr->second=DigitizeEnergy(itr->first->GetEnergy());
		Xtopsum+=itr->first->GetX()*itr->second;	
		Ytopsum+=itr->first->GetY()*itr->second;	
		energysum+=itr->second;	
		std::cout << "HITLAYER " << itr->first->GetLayer() << std::endl;
	}		
	barypositions={Xtopsum/energysum,Ytopsum/energysum,Zbary};
}



void splitcalHPLCluster::Clusterize(){
	//pick seed for Cluster
	splitcalHit *Seedhit=nullptr;
	StopSignal = kFALSE;
//	std::cout << "0 Clusterize" << std::endl;
	while(StopSignal==kFALSE){
		StopSignal = kTRUE;
//		std::cout << "1 Clusterize" << std::endl;
		for(int i=0;i<_vectorOfHits.size();i++){
			if(!(_vectorOfHits.at(i)->IsClustered()) && i==0){
//				std::cout << "2 Clusterize" << std::endl;
				Seedhit = _vectorOfHits.at(i);
				vectorOfClusters.push_back({_vectorOfHits.at(i),0});
				_vectorOfHits.at(i)->SetIsClustered(1);
				StopSignal=kFALSE;
				}
			if(!(_vectorOfHits.at(i)->IsClustered()) && TMath::Sqrt((Seedhit->GetX()-_vectorOfHits.at(i)->GetX())*(Seedhit->GetX()-_vectorOfHits.at(i)->GetX())+(Seedhit->GetY()-_vectorOfHits.at(i)->GetY())*(Seedhit->GetY()-_vectorOfHits.at(i)->GetY())) <10 && TMath::Sqrt(_vectorOfHits.at(i)->GetZ() - Seedhit->GetZ()) <1){
				vectorOfClusters.push_back({_vectorOfHits.at(i),0});
//				std::cout << "3 Clusterize" << std::endl;
				_vectorOfHits.at(i)->SetIsClustered(1);
//				std::cout << "4 Clusterize" << std::endl;
				StopSignal = kFALSE;
			}
		}
	}
//	std::cout << "5 Clusterize" << std::endl;
	Analyze();
//	std::cout << "ENDCLUS" << std::endl;
}



// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
splitcalHPLCluster::~splitcalHPLCluster() { }
// -------------------------------------------------------------------------

// -----   Public method Print   -------------------------------------------
void splitcalHPLCluster::Print() const
{

  std::cout<< "-I- splitcalHPLCluster: " <<std::endl;
  std::cout<< "    (eta,phi,energy) = " 
	   << _eta << " ,  "
	   << _phi << " ,  "
	   << _energy << std::endl;  
  std::cout<< "    start(x,y,z) = " 
	   << _start.X() << " ,  "
	   << _start.Y() << " ,  "
	   << _start.Z() << std::endl;  
  std::cout<< "    end(x,y,z) = " 
	   << _end.X() << " ,  "
	   << _end.Y() << " ,  "
	   << _end.Z() << std::endl;  
   std::cout<< "------- " <<std::endl;
}
// -------------------------------------------------------------------------

ClassImp(splitcalHPLCluster)

