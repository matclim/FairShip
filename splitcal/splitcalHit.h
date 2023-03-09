#ifndef SPLITCALHIT_H
#define SPLITCALHIT_H 1


#include "ShipHit.h"
#include "splitcalPoint.h"
#include "TObject.h"
#include "TVector3.h"

class splitcalHit : public ShipHit
{
  public:

    /** Default constructor **/
    splitcalHit();

    /** Constructor with arguments
     *@param detID    Detector ID
     *@param digi      digitized/measured TDC 
     *@param flag      True/False, false if there is another hit with smaller tdc 
     **/
    splitcalHit(Int_t detID, Float_t tdc);
    splitcalHit(splitcalPoint* p, Double_t t0);
    splitcalHit(double recoX, double recoY, double recoZ, double recoE, double pex, double pey);
    /** Destructor **/
    virtual ~splitcalHit();

    /** Output to screen **/
    virtual void Print() const;
    Float_t GetTDC() const {return fdigi;}
    void setInvalid() {flag = false;}
    bool isValid() const {return flag;}

    std::string GetPaddedString(int& id); 
    std::string GetDetectorElementName(int isPrec ,int& id);
    void Decoder(int& id, int& nLayer, int& nModuleX, int& nMdouleY, int& nStrip);
    void Decoder(std::string& encodedID, int& nLayer, int& nModuleX, int& nMdouleY, int& nStrip);

    void SetXYZ(double& x, double& y, double& z) {_x = x; _y = y; _z = z;}
    void SetIDs(int& isPrecision, int& nLayer, int& nModuleX, int& nModuleY, int& nStrip) {_isPrecisionLayer = isPrecision; _nLayer = nLayer; _nModuleX = nModuleX; _nModuleY=nModuleY; _nStrip = nStrip;}
    void SetisPrecisionLayer(int isPrecision){_isPrecisionLayer=isPrecision;} 
    void SetEnergy(double& e) {_energy = e;}
    void UpdateEnergy(double e) {_energy = _energy+e;}
    void SetIsX(bool x) {_isX = x;}
    void SetIsY(bool y) {_isY = y;}
    void SetIsUsed(int u) {_isUsed = u;}
    void SetIsClustered(int u) {_isClustered = u;}
    void SetXYZErrors(double xError, double yError, double zError) {_xError = xError; _yError = yError; _zError = zError;}
    /* void SetClusterIndex(int i) {_clusterIndex = i;} */
    /* void SetEnergyWeight(double w) {_energyWeight = w;} */
    void AddClusterIndex(int i) {_vecClusterIndices.push_back(i);}
    void AddEnergyWeight(double w) {_vecEnergyWeights.push_back(w);}

    double GetX() {return _x;}
    double GetY() {return _y;}
    double GetZ() {return _z;}
    double GetEnergy() {return _energy;}
    double GetEnergyForCluster(int i);
    int GetIsPrecisionLayer() {return _isPrecisionLayer;}
    int GetLayerNumber() {return _nLayer;}
    int GetModuleXNumber() {return _nModuleX;}
    int GetModuleYNumber() {return _nModuleY;}
    int GetStripNumber() {return _nStrip;}
    bool IsX() {return _isX;}
    bool IsY() {return _isY;}
    int IsUsed() {return _isUsed;}
    int IsClustered() {return _isClustered;}
    double GetXError() {return _xError;}
    double GetYError() {return _yError;}
    double GetZError() {return _zError;}
    /* double GetEnergyWeight() {return _energyWeight;} */
    /* int GetClusterIndex() {return _clusterIndex;} */
    std::vector<int > GetClusterIndices() {return _vecClusterIndices;}
    std::vector<double > GetEnergyWeights() {return _vecEnergyWeights;}
    bool IsShared() {return GetClusterIndices().size()>1; }
    double GetEnergyWeightForIndex(int index);
    int GetDetectorID(){return _detID;};
    /** Copy constructor **/
    splitcalHit(const splitcalHit& point);
    splitcalHit operator=(const splitcalHit& point);
    std::vector<splitcalHit*> GetvecofExtraHits(){return _vecofExtraHits;};
    Double_t GetEnergyFraction(){return _energyfraction;};

    /*for testing*/

    Bool_t firstextraflag=kFALSE; // for testing
    Int_t GetsizeofExtraHits(){return sizeofExtraHits;};
    Int_t GetsizeofExtraDiffusionHits(){return sizeofExtraDiffusionHits;};
    Double_t GetPx(){return px;};
    Double_t GetPy(){return py;};
    Int_t GetLayer(){return HPLayer;}
  private:

    /*for testing*/
    int sizeofExtraHits; 
    int sizeofExtraDiffusionHits=0; 
    double px,py;

    Float_t flag;   
    double _x, _y, _z, _xError, _yError, _zError;

    double _energy,_cloudsigmaX=0.00068,_cloudsigmaY=0.00068,HPLentresh=0.000001;// _energyWeight; /*cloudsigmas to be modified where needed, perhaps introduce an energy dependence ?*/
    /* std::string _nameSubDetector;  */
    int _isPrecisionLayer, _nLayer, _nModuleX, _nModuleY, _nStrip, _isUsed, _detID, _isClustered=0;// _clusterIndex;
    bool _isX, _isY;
    std::vector<double > _vecEnergyWeights;
    std::vector<int > _vecClusterIndices;
    std::vector<splitcalHit*> _vecofExtraHits;
    void ComputeHPL(double pointX,double pointY,double pointZ,double pointE,TVector3 *pointPVector,int _detID);
    void ComputeHPLInterim(double pointX,double pointY,double pointZ);
    void ComputeCloudSigma(double pointE);
    void SelectHPLayer(int DETID);
    std::array<double,2> ComputeBaryCenter(std::vector<std::array<double,4>> vecofrec,std::vector<double> vecofe);
    std::vector<double> DigitizeEnergy(std::vector<double> vecofe);
    Double_t ComputeHPLDiffusedEnergy(double x1, double x2, double y1, double y2, double mux, double muy);
    Double_t _energyfraction;
    int HPLayer=-1;
    ClassDef(splitcalHit,3);

};

#endif
