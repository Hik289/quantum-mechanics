/*Basic Markov chain Monte Carlo for nematicity model.*/

#ifndef _BSCMC_HPP
#define _BSCMC_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include "data_analyzer.hpp"
#include "Eigen/Eigen"
#include "model.hpp"
#include "random_Mersenne_twister.hpp"

#define PI 3.14159265358979323846

using std::vector;
using std::string;
using std::ofstream;
using std::ifstream;
using std::ios;
using Eigen::VectorXd;
using Eigen::Matrix;
using Eigen::SparseMatrix;
using Eigen::Dynamic;
using Eigen::MatrixXd;
using Eigen::Triplet;
using Eigen::ColPivHouseholderQR;
using Eigen::JacobiSVD;
using Eigen::ComputeThinU;
using Eigen::ComputeThinV;
using Eigen::Upper;

Matrix<double,4,4> genFourSiteMatExp(const double x,const double y,const double m);
void udvPivQR(const MatrixXd& mat,MatrixXd& matU,SparseMatrix<double>& matD,MatrixXd& matV);
void udvSinVal(const MatrixXd& mat,MatrixXd& matU,SparseMatrix<double>& matD,MatrixXd& matV);
void calcUdv(const MatrixXd& mat,MatrixXd& matU,SparseMatrix<double>& matD,MatrixXd& matV);

class MatK{
    public:

    SparseMatrix<double> A1PQuar;
    SparseMatrix<double> A2PQuar;
    SparseMatrix<double> A1PHalf;
    SparseMatrix<double> A2PHalf;
    SparseMatrix<double> A1P;
    SparseMatrix<double> A2P;

    SparseMatrix<double> A1NQuar;
    SparseMatrix<double> A2NQuar;
    SparseMatrix<double> A1NHalf;
    SparseMatrix<double> A2NHalf;
    SparseMatrix<double> A1N;
    SparseMatrix<double> A2N;

    SparseMatrix<double> B1PQuar;
    SparseMatrix<double> B2PQuar;
    SparseMatrix<double> B1PHalf;
    SparseMatrix<double> B2PHalf;
    SparseMatrix<double> B1P;
    SparseMatrix<double> B2P;

    SparseMatrix<double> B1NQuar;
    SparseMatrix<double> B2NQuar;
    SparseMatrix<double> B1NHalf;
    SparseMatrix<double> B2NHalf;
    SparseMatrix<double> B1N;
    SparseMatrix<double> B2N;
};

class GrnFunc{
    public:
    MatrixXd equA;
    MatrixXd equB;
    
    MatrixXd modEquA;
    MatrixXd modEquB;
};

class UdvRcrd{
    public:
    /*X-X-XX*/
    /*orbital (A or B)-udv decomposition components (U or D or V)-TauZero or BetaTau*/

    MatrixXd* AUTZ;
    SparseMatrix<double>* ADTZ;
    MatrixXd* AVTZ;

    MatrixXd* BUTZ;
    SparseMatrix<double>* BDTZ;
    MatrixXd* BVTZ;

    MatrixXd* AUBT;
    SparseMatrix<double>* ADBT;
    MatrixXd* AVBT;

    MatrixXd* BUBT;
    SparseMatrix<double>* BDBT;
    MatrixXd* BVBT;

    /*construct and destruct*/
    UdvRcrd(const int numUdvIn=1);
    UdvRcrd(const UdvRcrd& obj);
    UdvRcrd& operator= (const UdvRcrd& obj);
    ~UdvRcrd();
    void setNumUdv(const int numUdvIn);
    private:
    int szArr;
};

class BscMC{
public:
    double invTem;
    double win;
    int numUdv;
    int udvIntv;

    /*non-free*/
    int numSlc;
    double timeIntv;
    double** conf;

    Model model;
    RanMT ranMT;
    MatK matK;
    GrnFunc grnFunc;
    UdvRcrd udvRcrd;

    /*********************************constructor and destructor*********************************/
    BscMC();
    BscMC(const BscMC& obj);
    BscMC& operator= (const BscMC& obj);
    ~BscMC();

    /*******************************parameter settings and showing*******************************/
    void setInvTem(const double invTemIn);
    void setWin(const double winIn);
    void setSeed(const int mtSeedIn);
    void setUdv(const int numUdvIn,const int udvIntvIn);
    void setLattSize(const int lenXIn,const int lenYIn);
    void setHmltParm(double* parmIn);

    void showParm(ostream& iout=cout)const;
    void showConf(ostream& iout=cout)const;

    void exportStts(string loct);
    void importStts(string loct);

    /***********************************configuration updating***********************************/
    void ranConf(const double ranMin=-10,const double ranMax=10);

    SparseMatrix<double> genMatExpK(const double termHopX,const double termHopY,const double termChem,const int CBInd);
    SparseMatrix<double> genMatExpV(const int slcInd,const int orbInd);
    void initMatK();

    void syncUdvRcrd();
    void syncGrnEquA(const int udvInd);
    void syncGrnEquB(const int udvInd);
    void syncGrnEqu(const int udvInd);
    void syncGrnModEquA(const int udvInd);
    void syncGrnModEquB(const int udvInd);
    void syncGrnModEqu(const int udvInd);

    void prepScanFromZero();
    void prepScanFromBeta();

    double wgtRatFer(const int siteUpd,const double stateChng)const;
    double wgtRatBos(const int timeUpd,const int siteUpd,const double stateChng)const;
    double wgtRat(const int timeUpd,const int siteUpd,const double stateChng)const;
    void updGrnModEqu(const int siteUpd,const double stateChng);

    void swpSlc(const int slcUpd);
    void swpSlc(const int slcUpd,double& accpRate);

    void evoGrnFwd(const int slcInd);
    void evoGrnBwd(const int slcInd);

    void scanZeroBeta();
    void scanBetaZero();
    void scanZeroBeta(double& accpRate,double& errGrn,double& condNum);
    void scanBetaZero(double& accpRate,double& errGrn,double& condNum);
    void cycScan(const int dir=0);
    void cycScan(double& accpRate,double& errGrn,double& condNum,const int dir=0);

    void printRoundOffInfo(ostream& iout=cout);
    double calcAccpRate();
    void adjWin();
    void detWin(const double winMin,const double winMax,ostream& iout=cout);

    /***********************************observable measurement***********************************/
    void mesrEquGrn(MatrixXd* resA,MatrixXd* resB);
    void mesrNeqGrn(MatrixXd* resA,MatrixXd* resB);

    void findSpacRSite(int rInd,int* res)const;
    void calcEquCorrSngSlc(MatrixXd& equGrnAIn,MatrixXd& equGrnBIn,vector<double>& res)const;
    void calcEquCorr(MatrixXd* equGrnAIn,MatrixXd* equGrnBIn,vector<double>& res)const;
    double mesrConjFld()const;
    double mesrConjFldSq()const;
    double mesrConjFldQt()const;
    double calcEquCorrFMSngSlc(MatrixXd& equGrnAIn,MatrixXd& equGrnBIn)const;
    double calcEquCorrFM(MatrixXd* equGrnAIn,MatrixXd* equGrnBIn)const;

    double calcOrdFMSngSlc(MatrixXd& equGrnAIn,MatrixXd& equGrnBIn)const;
    double calcOrdFM(MatrixXd* equGrnAIn,MatrixXd* equGrnBIn)const;
    double calcOrdFMSq(MatrixXd* equGrnAIn,MatrixXd* equGrnBIn)const;
    double calcOrdFMQt(MatrixXd* equGrnAIn,MatrixXd* equGrnBIn)const;
    double calcPairCorrS(MatrixXd* equGrnAIn,MatrixXd* equGrnBIn)const;

    double calcAvrgOcupSngSlc(MatrixXd& equGrnIn)const;
    double calcAvrgOcup(MatrixXd* equGrnIn)const;

    void calcNeqGrnKSngOrbt(MatrixXd* neqOrigR,double** neqK)const;
    void calcNeqGrnK(MatrixXd* neqOrigRA,MatrixXd* neqOrigRB,double** neqKA,double** neqKB)const;
    /*calculate pseudo density of states*/
    double calcPDOS(double** neqKIn)const;

    /*calculate and show occupation configuration*/
    void showOcupConf(const string filePath);
};

/*used for checkerboard decomposition*/
/*generate exponential form of four-site matrix*/
/*input: hopping along x and y (x,y) as well as chemical potential (m)*/
/*output: 4 by 4 double-type-matrix*/
Matrix<double,4,4> genFourSiteMatExp(const double x,const double y,const double m){
    Matrix<double,4,4> res;
    double chX=cosh(x);
    double chY=cosh(y);
    double chMH=cosh(m/2);
    double shX=sinh(x);
    double shY=sinh(y);
    double shMH=sinh(m/2);

    res<<chX*chY*(chMH-shMH),shX*chY*(chMH-shMH),shX*shY*(chMH-shMH),chX*shY*(chMH-shMH),
         shX*chY*(chMH-shMH),chX*chY*(chMH-shMH),chX*shY*(chMH-shMH),shX*shY*(chMH-shMH),
         shX*shY*(chMH-shMH),chX*shY*(chMH-shMH),chX*chY*(chMH-shMH),shX*chY*(chMH-shMH),
         chX*shY*(chMH-shMH),shX*shY*(chMH-shMH),shX*chY*(chMH-shMH),chX*chY*(chMH-shMH);
    return res;
}

/*UDV decomposition using column-pivoting QR algorithm*/
void udvPivQR(const MatrixXd& mat,MatrixXd& matU,SparseMatrix<double>& matD,MatrixXd& matV){
    ColPivHouseholderQR< MatrixXd > pivQR;
    pivQR.compute(mat);
    MatrixXd matP=pivQR.colsPermutation();
    MatrixXd matR=pivQR.matrixR().triangularView<Upper>();

    matU=pivQR.matrixQ();
    const int dim=mat.cols();
    vector< Triplet<double> > triList(dim);
    for(int i=0;i<dim;i++){
        triList[i]=Triplet<double>(i,i,matR(i,i));
    }
    matD.setFromTriplets(triList.begin(),triList.end());
    for(int i=0;i<dim;i++){
        triList[i]=Triplet<double>(i,i,1/matR(i,i));
    }
    SparseMatrix<double> matDInv(dim,dim);
    matDInv.setFromTriplets(triList.begin(),triList.end());
    matV=matR*matP.transpose();
    matV=matDInv*matV;
}

/*UDV decomposition using two-sided Jacobi transformation*/
void udvSinVal(const MatrixXd& mat,MatrixXd& matU,SparseMatrix<double>& matD,MatrixXd& matV){
    const int dim=mat.cols();
    JacobiSVD< MatrixXd > JacSVD;
    JacSVD.compute(mat,ComputeThinU | ComputeThinV);
    matU=JacSVD.matrixU();
    matV=JacSVD.matrixV().transpose();
    MatrixXd sinValList=JacSVD.singularValues();
    vector< Triplet<double> >triList(dim);
    for(int i=0;i<dim;i++){
        triList[i]=Triplet<double>(i,i,sinValList(i));
    }
    matD.setFromTriplets(triList.begin(),triList.end());
}

/*choose which method to use for UDV factorization*/
/*for large matrice pivoting QR is slower but more accurate than Jacobi SVD*/
void calcUdv(const MatrixXd& mat,MatrixXd& matU,SparseMatrix<double>& matD,MatrixXd& matV){
    udvPivQR(mat,matU,matD,matV);
    // udvSinVal(mat,matU,matD,matV);
}

UdvRcrd::UdvRcrd(const int numUdvIn){
    this->szArr=numUdvIn;
    this->AUTZ=new MatrixXd[numUdvIn];
    this->BUTZ=new MatrixXd[numUdvIn];
    this->AUBT=new MatrixXd[numUdvIn];
    this->BUBT=new MatrixXd[numUdvIn];

    this->ADTZ=new SparseMatrix<double>[numUdvIn];
    this->BDTZ=new SparseMatrix<double>[numUdvIn];
    this->ADBT=new SparseMatrix<double>[numUdvIn];
    this->BDBT=new SparseMatrix<double>[numUdvIn];

    this->AVTZ=new MatrixXd[numUdvIn];
    this->BVTZ=new MatrixXd[numUdvIn];
    this->AVBT=new MatrixXd[numUdvIn];
    this->BVBT=new MatrixXd[numUdvIn];
}

UdvRcrd::~UdvRcrd(){
    delete[] this->AUTZ;
    delete[] this->BUTZ;
    delete[] this->AUBT;
    delete[] this->BUBT;

    delete[] this->ADTZ;
    delete[] this->BDTZ;
    delete[] this->ADBT;
    delete[] this->BDBT;

    delete[] this->AVTZ;
    delete[] this->BVTZ;
    delete[] this->AVBT;
    delete[] this->BVBT;
}

UdvRcrd::UdvRcrd(const UdvRcrd& obj){
    this->szArr=obj.szArr;
    this->AUTZ=new MatrixXd[this->szArr];
    this->BUTZ=new MatrixXd[this->szArr];
    this->AUBT=new MatrixXd[this->szArr];
    this->BUBT=new MatrixXd[this->szArr];

    this->ADTZ=new SparseMatrix<double>[this->szArr];
    this->BDTZ=new SparseMatrix<double>[this->szArr];
    this->ADBT=new SparseMatrix<double>[this->szArr];
    this->BDBT=new SparseMatrix<double>[this->szArr];

    this->AVTZ=new MatrixXd[this->szArr];
    this->BVTZ=new MatrixXd[this->szArr];
    this->AVBT=new MatrixXd[this->szArr];
    this->BVBT=new MatrixXd[this->szArr];

    for(int i=0;i<this->szArr;i++){
        this->AUTZ[i]=obj.AUTZ[i];
        this->BUTZ[i]=obj.BUTZ[i];
        this->AUBT[i]=obj.AUBT[i];
        this->BUBT[i]=obj.BUBT[i];

        this->ADTZ[i]=obj.ADTZ[i];
        this->BDTZ[i]=obj.BDTZ[i];
        this->ADBT[i]=obj.ADBT[i];
        this->BDBT[i]=obj.BDBT[i];

        this->AVTZ[i]=obj.AVTZ[i];
        this->BVTZ[i]=obj.BVTZ[i];
        this->AVBT[i]=obj.AVBT[i];
        this->BVBT[i]=obj.BVBT[i];
    }
}

UdvRcrd& UdvRcrd::operator= (const UdvRcrd& obj){
    if(this==&obj){
        return *this;
    }
    delete[] this->AUTZ;
    delete[] this->BUTZ;
    delete[] this->AUBT;
    delete[] this->BUBT;

    delete[] this->ADTZ;
    delete[] this->BDTZ;
    delete[] this->ADBT;
    delete[] this->BDBT;

    delete[] this->AVTZ;
    delete[] this->BVTZ;
    delete[] this->AVBT;
    delete[] this->BVBT;

    this->szArr=obj.szArr;
    this->AUTZ=new MatrixXd[this->szArr];
    this->BUTZ=new MatrixXd[this->szArr];
    this->AUBT=new MatrixXd[this->szArr];
    this->BUBT=new MatrixXd[this->szArr];

    this->ADTZ=new SparseMatrix<double>[this->szArr];
    this->BDTZ=new SparseMatrix<double>[this->szArr];
    this->ADBT=new SparseMatrix<double>[this->szArr];
    this->BDBT=new SparseMatrix<double>[this->szArr];

    this->AVTZ=new MatrixXd[this->szArr];
    this->BVTZ=new MatrixXd[this->szArr];
    this->AVBT=new MatrixXd[this->szArr];
    this->BVBT=new MatrixXd[this->szArr];

    for(int i=0;i<this->szArr;i++){
        this->AUTZ[i]=obj.AUTZ[i];
        this->BUTZ[i]=obj.BUTZ[i];
        this->AUBT[i]=obj.AUBT[i];
        this->BUBT[i]=obj.BUBT[i];

        this->ADTZ[i]=obj.ADTZ[i];
        this->BDTZ[i]=obj.BDTZ[i];
        this->ADBT[i]=obj.ADBT[i];
        this->BDBT[i]=obj.BDBT[i];

        this->AVTZ[i]=obj.AVTZ[i];
        this->BVTZ[i]=obj.BVTZ[i];
        this->AVBT[i]=obj.AVBT[i];
        this->BVBT[i]=obj.BVBT[i];
    }
    return *this;
}

void UdvRcrd::setNumUdv(const int numUdvIn){
    delete[] this->AUTZ;
    delete[] this->BUTZ;
    delete[] this->AUBT;
    delete[] this->BUBT;

    delete[] this->ADTZ;
    delete[] this->BDTZ;
    delete[] this->ADBT;
    delete[] this->BDBT;

    delete[] this->AVTZ;
    delete[] this->BVTZ;
    delete[] this->AVBT;
    delete[] this->BVBT;

    this->szArr=numUdvIn;
    this->AUTZ=new MatrixXd[this->szArr];
    this->BUTZ=new MatrixXd[this->szArr];
    this->AUBT=new MatrixXd[this->szArr];
    this->BUBT=new MatrixXd[this->szArr];

    this->ADTZ=new SparseMatrix<double>[this->szArr];
    this->BDTZ=new SparseMatrix<double>[this->szArr];
    this->ADBT=new SparseMatrix<double>[this->szArr];
    this->BDBT=new SparseMatrix<double>[this->szArr];

    this->AVTZ=new MatrixXd[this->szArr];
    this->BVTZ=new MatrixXd[this->szArr];
    this->AVBT=new MatrixXd[this->szArr];
    this->BVBT=new MatrixXd[this->szArr];
}

BscMC::BscMC(){
    this->invTem=1;
    this->win=1;
    this->ranMT.setSeed(0);
    this->numUdv=1;
    this->udvIntv=1;
    this->numSlc=this->numUdv*this->udvIntv;
    this->timeIntv=this->invTem/this->numSlc;
    this->model.setLattSize(4,4);
    double tmp[NUM_HMLT_PARM];
    for(int i=0;i<NUM_HMLT_PARM;i++){
        tmp[i]=0;
    }
    this->model.setHmltParm(tmp);
    udvRcrd.setNumUdv(this->numUdv);

    this->conf=new double*[this->numSlc];
    for(int i=0;i<this->numSlc;i++){
        this->conf[i]=new double[this->model.numSite];
        for(int j=0;j<this->model.numSite;j++){
            this->conf[i][j]=0;
        }
    }
}

BscMC::~BscMC(){
    for(int i=0;i<this->numSlc;i++){
        delete[] this->conf[i];
    }
    delete[] this->conf;
}

BscMC::BscMC(const BscMC& obj){
    this->invTem=obj.invTem;
    this->win=obj.win;
    this->numUdv=obj.numUdv;
    this->udvIntv=obj.udvIntv;
    this->numSlc=this->numUdv*this->udvIntv;
    this->timeIntv=this->invTem/this->numSlc;
    this->model=obj.model;
    this->ranMT=obj.ranMT;
    this->matK=obj.matK;
    this->grnFunc=obj.grnFunc;
    this->udvRcrd=obj.udvRcrd;

    this->conf=new double*[this->numSlc];
    for(int i=0;i<this->numSlc;i++){
        this->conf[i]=new double[this->model.numSite];
        for(int j=0;j<this->model.numSite;j++){
            this->conf[i][j]=obj.conf[i][j];
        }
    }
}

BscMC& BscMC::operator= (const BscMC& obj){
    if(this==&obj){
        return *this;
    }
    for(int i=0;i<this->numSlc;i++){
        delete[] this->conf[i];
    }
    delete[] this->conf;

    this->invTem=obj.invTem;
    this->win=obj.win;
    this->numUdv=obj.numUdv;
    this->udvIntv=obj.udvIntv;
    this->numSlc=this->numUdv*this->udvIntv;
    this->timeIntv=this->invTem/this->numSlc;
    this->model=obj.model;
    this->ranMT=obj.ranMT;
    this->matK=obj.matK;
    this->grnFunc=obj.grnFunc;
    this->udvRcrd=obj.udvRcrd;

    this->conf=new double*[this->numSlc];
    for(int i=0;i<this->numSlc;i++){
        this->conf[i]=new double[this->model.numSite];
        for(int j=0;j<this->model.numSite;j++){
            this->conf[i][j]=obj.conf[i][j];
        }
    }

    return *this;
}

void BscMC::setInvTem(const double invTemIn){
    this->invTem=invTemIn;
    this->timeIntv=this->invTem/this->numSlc;
}
void BscMC::setWin(const double winIn){
    this->win=winIn;
}
void BscMC::setSeed(const int mtSeedIn){
    this->ranMT.setSeed(mtSeedIn);
}
void BscMC::setUdv(const int numUdvIn,const int udvIntvIn){
    for(int i=0;i<this->numSlc;i++){
        delete[] this->conf[i];
    }
    delete[] this->conf;

    this->numUdv=numUdvIn;
    this->udvIntv=udvIntvIn;
    this->numSlc=this->numUdv*this->udvIntv;
    this->timeIntv=this->invTem/this->numSlc;
    this->udvRcrd.setNumUdv(this->numUdv);

    this->conf=new double*[this->numSlc];
    for(int i=0;i<this->numSlc;i++){
        this->conf[i]=new double[this->model.numSite];
        for(int j=0;j<this->model.numSite;j++){
            this->conf[i][j]=0;
        }
    }
}

void BscMC::setLattSize(const int lenXIn,const int lenYIn){
    for(int i=0;i<this->numSlc;i++){
        delete[] this->conf[i];
    }
    delete[] this->conf;

    this->model.setLattSize(lenXIn,lenYIn);
    this->conf=new double*[this->numSlc];

    for(int i=0;i<this->numSlc;i++){
        this->conf[i]=new double[this->model.numSite];
        for(int j=0;j<this->model.numSite;j++){
            this->conf[i][j]=0;
        }
    }
}

void BscMC::setHmltParm(double* parmIn){
    this->model.setHmltParm(parmIn);
}

void BscMC::showParm(ostream& iout)const{
    iout<<"===========Monte Carlo paramters=============="<<endl;
    iout<<"lattice size: "<<this->model.lenX<<"x"<<this->model.lenY<<endl;
    iout<<"inverse Tmperature: "<<this->invTem<<endl;
    iout<<"number of time slice: "<<this->numSlc<<endl;
    iout<<"number of udv decomposition: "<<this->numUdv<<endl;
    iout<<"interval of udv decomposition: "<<this->udvIntv<<endl;
    iout<<"time interval: "<<this->timeIntv<<endl;
    iout<<"random seed: "<<this->ranMT.readSeed()<<endl;
    iout<<"hopping strength:        "<<this->model.strHop<<endl;
    iout<<"hopping ratio:           "<<this->model.ratHop<<endl;
    iout<<"coupling strength:       "<<this->model.strCpl<<endl;
    iout<<"dynamic strength:        "<<this->model.strDyn<<endl;
    iout<<"gradient strength:       "<<this->model.strGrad<<endl;
    iout<<"squared term strength:   "<<this->model.strSq<<endl;
    iout<<"quartic term strength:   "<<this->model.strQt<<endl;
    iout<<"chemical potential:      "<<this->model.chem<<endl;
    iout<<"chemical potential bias: "<<this->model.biasChem<<endl;
}
void BscMC::showConf(ostream& iout)const{
    iout<<"===========current configuration=============="<<endl;
    for(int slcInd=0;slcInd<this->numSlc;slcInd++){
        iout<<"time slice "<<slcInd+1<<endl;
        for(int y=0;y<this->model.lenY;y++){
            for(int x=0;x<this->model.lenX;x++){
                iout<<this->conf[slcInd][y*this->model.lenX+x]<<" ";
            }
            iout<<endl;
        }
    }
}

void BscMC::ranConf(const double ranMin,const double ranMax){
    for(int slcInd=0;slcInd<this->numSlc;slcInd++){
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            this->conf[slcInd][siteInd]=this->ranMT.ranReal(ranMin,ranMax);
        }
    }
}

/*used to generate K-matrice*/
/*time interval is absorbed into hoppings*/
/*termChem: chemical potential*/
/*CBInd means checkerboard index (1 or 2)*/
SparseMatrix<double> BscMC::genMatExpK(const double termHopX,const double termHopY,const double termChem,const int CBInd){
    int** FSList;
    if(CBInd==1){
        FSList=this->model.fourSite1;
    }
    if(CBInd==2){
        FSList=this->model.fourSite2;
    }

    Matrix<double,4,4> fourSiteMatExp=genFourSiteMatExp(termHopX,termHopY,termChem);
    vector< Triplet<double> > triList(4*this->model.numSite);
    for(int FSInd=0;FSInd<(this->model.numSite/4);FSInd++){
        for(int i=0;i<4;i++){
            for(int j=0;j<4;j++){
                int triInd=FSInd*16+i*4+j;
                triList[triInd]=Triplet<double>(FSList[FSInd][i],FSList[FSInd][j],fourSiteMatExp(i,j));
            }
        }
    }
    SparseMatrix<double> res(this->model.numSite,this->model.numSite);
    res.setFromTriplets(triList.begin(),triList.end());
    return res;
}

/*basic matrice generation*/
void BscMC::initMatK(){
    double termHopX,termHopY,termChem;

    termHopX=(0.25*this->timeIntv)*(-this->model.strHop);
    termHopY=(0.25*this->timeIntv)*(this->model.ratHop*this->model.strHop);
    termChem=(0.25*this->timeIntv)*(this->model.chem+this->model.biasChem);
    this->matK.A1PQuar=genMatExpK(termHopX,termHopY,termChem,1);
    this->matK.A2PQuar=genMatExpK(termHopX,termHopY,termChem,2);

    termHopX=(0.5*this->timeIntv)*(-this->model.strHop);
    termHopY=(0.5*this->timeIntv)*(this->model.ratHop*this->model.strHop);
    termChem=(0.5*this->timeIntv)*(this->model.chem+this->model.biasChem);
    this->matK.A1PHalf=genMatExpK(termHopX,termHopY,termChem,1);
    this->matK.A2PHalf=genMatExpK(termHopX,termHopY,termChem,2);

    termHopX=(1.0*this->timeIntv)*(-this->model.strHop);
    termHopY=(1.0*this->timeIntv)*(this->model.ratHop*this->model.strHop);
    termChem=(1.0*this->timeIntv)*(this->model.chem+this->model.biasChem);
    this->matK.A1P=genMatExpK(termHopX,termHopY,termChem,1);
    this->matK.A2P=genMatExpK(termHopX,termHopY,termChem,2);
    /************************************/
    termHopX=(-0.25*this->timeIntv)*(-this->model.strHop);
    termHopY=(-0.25*this->timeIntv)*(this->model.ratHop*this->model.strHop);
    termChem=(-0.25*this->timeIntv)*(this->model.chem+this->model.biasChem);
    this->matK.A1NQuar=genMatExpK(termHopX,termHopY,termChem,1);
    this->matK.A2NQuar=genMatExpK(termHopX,termHopY,termChem,2);

    termHopX=(-0.5*this->timeIntv)*(-this->model.strHop);
    termHopY=(-0.5*this->timeIntv)*(this->model.ratHop*this->model.strHop);
    termChem=(-0.5*this->timeIntv)*(this->model.chem+this->model.biasChem);
    this->matK.A1NHalf=genMatExpK(termHopX,termHopY,termChem,1);
    this->matK.A2NHalf=genMatExpK(termHopX,termHopY,termChem,2);

    termHopX=(-1.0*this->timeIntv)*(-this->model.strHop);
    termHopY=(-1.0*this->timeIntv)*(this->model.ratHop*this->model.strHop);
    termChem=(-1.0*this->timeIntv)*(this->model.chem+this->model.biasChem);
    this->matK.A1N=genMatExpK(termHopX,termHopY,termChem,1);
    this->matK.A2N=genMatExpK(termHopX,termHopY,termChem,2);
    // /************************************/
    termHopX=(0.25*this->timeIntv)*(this->model.ratHop*this->model.strHop);
    termHopY=(0.25*this->timeIntv)*(-this->model.strHop);
    termChem=(0.25*this->timeIntv)*(this->model.chem-this->model.biasChem);
    this->matK.B1PQuar=genMatExpK(termHopX,termHopY,termChem,1);
    this->matK.B2PQuar=genMatExpK(termHopX,termHopY,termChem,2);

    termHopX=(0.5*this->timeIntv)*(this->model.ratHop*this->model.strHop);
    termHopY=(0.5*this->timeIntv)*(-this->model.strHop);
    termChem=(0.5*this->timeIntv)*(this->model.chem-this->model.biasChem);
    this->matK.B1PHalf=genMatExpK(termHopX,termHopY,termChem,1);
    this->matK.B2PHalf=genMatExpK(termHopX,termHopY,termChem,2);

    termHopX=(1.0*this->timeIntv)*(this->model.ratHop*this->model.strHop);
    termHopY=(1.0*this->timeIntv)*(-this->model.strHop);
    termChem=(1.0*this->timeIntv)*(this->model.chem-this->model.biasChem);
    this->matK.B1P=genMatExpK(termHopX,termHopY,termChem,1);
    this->matK.B2P=genMatExpK(termHopX,termHopY,termChem,2);
    // /************************************/
    termHopX=(-0.25*this->timeIntv)*(this->model.ratHop*this->model.strHop);
    termHopY=(-0.25*this->timeIntv)*(-this->model.strHop);
    termChem=(-0.25*this->timeIntv)*(this->model.chem-this->model.biasChem);
    this->matK.B1NQuar=genMatExpK(termHopX,termHopY,termChem,1);
    this->matK.B2NQuar=genMatExpK(termHopX,termHopY,termChem,2);

    termHopX=(-0.5*this->timeIntv)*(this->model.ratHop*this->model.strHop);
    termHopY=(-0.5*this->timeIntv)*(-this->model.strHop);
    termChem=(-0.5*this->timeIntv)*(this->model.chem-this->model.biasChem);
    this->matK.B1NHalf=genMatExpK(termHopX,termHopY,termChem,1);
    this->matK.B2NHalf=genMatExpK(termHopX,termHopY,termChem,2);

    termHopX=(-1.0*this->timeIntv)*(this->model.ratHop*this->model.strHop);
    termHopY=(-1.0*this->timeIntv)*(-this->model.strHop);
    termChem=(-1.0*this->timeIntv)*(this->model.chem-this->model.biasChem);
    this->matK.B1N=genMatExpK(termHopX,termHopY,termChem,1);
    this->matK.B2N=genMatExpK(termHopX,termHopY,termChem,2);
}

/*given configuration and time slice, generate V-matrice*/
SparseMatrix<double> BscMC::genMatExpV(const int slcInd,const int orbInd){
    double preFac;
    if(orbInd==0){
        preFac=-this->timeIntv*this->model.strCpl;
    }else{
        preFac=this->timeIntv*this->model.strCpl;
    }
    vector< Triplet<double> > triList(this->model.numSite);
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        triList[siteInd]=Triplet<double>(siteInd,siteInd,exp(preFac*this->conf[slcInd][siteInd]));
    }
    SparseMatrix<double> res(this->model.numSite,this->model.numSite);
    res.setFromTriplets(triList.begin(),triList.end());
    return res;
}

/*based on configuration, calculate B-matrice stored as UDV-form*/
void BscMC::syncUdvRcrd(){
    /*containers initialize*/
    MatrixXd matU(this->model.numSite,this->model.numSite);
    MatrixXd matV(this->model.numSite,this->model.numSite);
    MatrixXd prod(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matD(this->model.numSite,this->model.numSite);

    MatrixXd matUTmp(this->model.numSite,this->model.numSite);
    MatrixXd matVTmp(this->model.numSite,this->model.numSite);
    MatrixXd matTmp(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDTmp(this->model.numSite,this->model.numSite);

    /*\tau-0 list for orbital A*/
    matU.setIdentity();
    matV.setIdentity();
    prod.setIdentity();
    matD.setIdentity();
    for(int slcInd=0;slcInd<this->numSlc;slcInd++){
        prod=this->matK.A1NQuar*prod;
        prod=this->matK.A2NHalf*prod;
        prod=this->matK.A1NQuar*prod;
        prod=this->genMatExpV(slcInd,0)*prod;
        prod=this->matK.A1NQuar*prod;
        prod=this->matK.A2NHalf*prod;
        prod=this->matK.A1NQuar*prod;
        if((slcInd+1)%this->udvIntv==0){
            matTmp=prod*matU;
            matTmp=matTmp*matD;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matU=matUTmp;
            matD=matDTmp;
            matV=matVTmp*matV;
            const int udvInd=(slcInd+1)/this->udvIntv-1;

            this->udvRcrd.AUTZ[udvInd]=matU;
            this->udvRcrd.ADTZ[udvInd]=matD;
            this->udvRcrd.AVTZ[udvInd]=matV;
            prod.setIdentity();
        }
    }

    /*\tau-0 list for orbital B*/
    matU.setIdentity();
    matV.setIdentity();
    prod.setIdentity();
    matD.setIdentity();
    for(int slcInd=0;slcInd<this->numSlc;slcInd++){
        prod=this->matK.B1NQuar*prod;
        prod=this->matK.B2NHalf*prod;
        prod=this->matK.B1NQuar*prod;
        prod=this->genMatExpV(slcInd,1)*prod;
        prod=this->matK.B1NQuar*prod;
        prod=this->matK.B2NHalf*prod;
        prod=this->matK.B1NQuar*prod;
        
        if((slcInd+1)%this->udvIntv==0){
            matTmp=prod*matU;
            matTmp=matTmp*matD;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matU=matUTmp;
            matD=matDTmp;
            matV=matVTmp*matV;
            const int udvInd=(slcInd+1)/this->udvIntv-1;

            this->udvRcrd.BUTZ[udvInd]=matU;
            this->udvRcrd.BDTZ[udvInd]=matD;
            this->udvRcrd.BVTZ[udvInd]=matV;
            prod.setIdentity();
        }
    }

    /*\beta-\tau list for orbital A*/
    matU.setIdentity();
    matV.setIdentity();
    prod.setIdentity();
    matD.setIdentity();
    for(int slcInd=this->numSlc-1;slcInd>=0;slcInd--){
        prod=this->matK.A1NQuar*prod;
        prod=this->matK.A2NHalf*prod;
        prod=this->matK.A1NQuar*prod;
        prod=this->genMatExpV(slcInd,0)*prod;
        prod=this->matK.A1NQuar*prod;
        prod=this->matK.A2NHalf*prod;
        prod=this->matK.A1NQuar*prod;
        
        if((this->numSlc-slcInd)%this->udvIntv==0){
            matTmp=prod*matU;
            matTmp=matTmp*matD;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matU=matUTmp;
            matD=matDTmp;
            matV=matVTmp*matV;
            const int udvInd=(this->numSlc-slcInd)/this->udvIntv-1;

            this->udvRcrd.AUBT[udvInd]=matU.transpose();
            this->udvRcrd.ADBT[udvInd]=matD;
            this->udvRcrd.AVBT[udvInd]=matV.transpose();
            prod.setIdentity();
        }
    }

    /*\beta-\tau list for orbital B*/
    matU.setIdentity();
    matV.setIdentity();
    prod.setIdentity();
    matD.setIdentity();
    for(int slcInd=this->numSlc-1;slcInd>=0;slcInd--){
        prod=this->matK.B1NQuar*prod;
        prod=this->matK.B2NHalf*prod;
        prod=this->matK.B1NQuar*prod;
        prod=this->genMatExpV(slcInd,1)*prod;
        prod=this->matK.B1NQuar*prod;
        prod=this->matK.B2NHalf*prod;
        prod=this->matK.B1NQuar*prod;
        
        if((this->numSlc-slcInd)%this->udvIntv==0){
            matTmp=prod*matU;
            matTmp=matTmp*matD;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matU=matUTmp;
            matD=matDTmp;
            matV=matVTmp*matV;
            const int udvInd=(this->numSlc-slcInd)/this->udvIntv-1;

            this->udvRcrd.BUBT[udvInd]=matU.transpose();
            this->udvRcrd.BDBT[udvInd]=matD;
            this->udvRcrd.BVBT[udvInd]=matV.transpose();
            prod.setIdentity();
        }
    }
}

/*calculate equal-time Green's function for a particular UDV index*/
/*based on B-matrice calculated by SysMC.calcAllUDVList()*/
/*UDV index can be -1,0,1,2,...,this->numUdv-1*/
/*for udvInd==-1, purely use data stored in matXXXBetaTauList*/
/*for udvInd==this->numUdv-1, purely use data stored in matXXXTauList*/
/*calculate Green's function for orbital A*/
void BscMC::syncGrnEquA(const int udvInd){
    vector< Triplet<double> > triList(this->model.numSite);
    MatrixXd matUR(this->model.numSite,this->model.numSite);
    MatrixXd matVR(this->model.numSite,this->model.numSite);
    vector<double> diagDR(this->model.numSite);
    vector<double> diagDRLargeInv(this->model.numSite);
    vector<double> diagDRSmall(this->model.numSite);
    SparseMatrix<double> matDRLargeInv(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDRSmall(this->model.numSite,this->model.numSite);
    if(udvInd==-1){
        matUR.setIdentity();
        matVR.setIdentity();
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            diagDR[siteInd]=1;
        }
    }else{
        matUR=this->udvRcrd.AUTZ[udvInd];
        matVR=this->udvRcrd.AVTZ[udvInd];
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            diagDR[siteInd]=this->udvRcrd.ADTZ[udvInd].coeffRef(siteInd,siteInd);
        }
    }
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        if(abs(diagDR[siteInd])>=1){
            diagDRLargeInv[siteInd]=1/diagDR[siteInd];
            diagDRSmall[siteInd]=1;
        }else{
            diagDRLargeInv[siteInd]=1;
            diagDRSmall[siteInd]=diagDR[siteInd];
        }
    }
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDRLargeInv[siteInd]);
    }
    matDRLargeInv.setFromTriplets(triList.begin(),triList.end());
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDRSmall[siteInd]);
    }
    matDRSmall.setFromTriplets(triList.begin(),triList.end());

    MatrixXd matUL(this->model.numSite,this->model.numSite);
    MatrixXd matVL(this->model.numSite,this->model.numSite);
    vector<double> diagDL(this->model.numSite);
    vector<double> diagDLLargeInv(this->model.numSite);
    vector<double> diagDLSmall(this->model.numSite);
    SparseMatrix<double> matDLLargeInv(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDLSmall(this->model.numSite,this->model.numSite);
    if(udvInd==this->numUdv-1){
        matUL.setIdentity();
        matVL.setIdentity();
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            diagDL[siteInd]=1;
        }
    }else{
        matUL=this->udvRcrd.AUBT[this->numUdv-2-udvInd];
        matVL=this->udvRcrd.AVBT[this->numUdv-2-udvInd];
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            diagDL[siteInd]=this->udvRcrd.ADBT[this->numUdv-2-udvInd].coeffRef(siteInd,siteInd);
        }
    }
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        if(abs(diagDL[siteInd])>=1){
            diagDLLargeInv[siteInd]=1/diagDL[siteInd];
            diagDLSmall[siteInd]=1;
        }else{
            diagDLLargeInv[siteInd]=1;
            diagDLSmall[siteInd]=diagDL[siteInd];
        }
    }
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDLLargeInv[siteInd]);
    }
    matDLLargeInv.setFromTriplets(triList.begin(),triList.end());
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDLSmall[siteInd]);
    }
    matDLSmall.setFromTriplets(triList.begin(),triList.end());

    MatrixXd matULInv=matUL.inverse();
    MatrixXd matURInv=matUR.inverse();
    MatrixXd matMid1,matMid2,matMid;
    matMid1=matURInv*matULInv;
    matMid1=matDRLargeInv*matMid1;
    matMid1=matMid1*matDLLargeInv;
    matMid2=matVR*matVL;
    matMid2=matDRSmall*matMid2;
    matMid2=matMid2*matDLSmall;
    matMid=matMid1+matMid2;
    matMid=matMid.inverse();
    this->grnFunc.equA=matMid;
    this->grnFunc.equA=matDLLargeInv*this->grnFunc.equA;
    this->grnFunc.equA=matULInv*this->grnFunc.equA;
    this->grnFunc.equA=this->grnFunc.equA*matDRLargeInv;
    this->grnFunc.equA=this->grnFunc.equA*matURInv;
}
void BscMC::syncGrnEquB(const int udvInd){
    vector< Triplet<double> > triList(this->model.numSite);
    MatrixXd matUR(this->model.numSite,this->model.numSite);
    MatrixXd matVR(this->model.numSite,this->model.numSite);
    vector<double> diagDR(this->model.numSite);
    vector<double> diagDRLargeInv(this->model.numSite);
    vector<double> diagDRSmall(this->model.numSite);
    SparseMatrix<double> matDRLargeInv(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDRSmall(this->model.numSite,this->model.numSite);
    if(udvInd==-1){
        matUR.setIdentity();
        matVR.setIdentity();
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            diagDR[siteInd]=1;
        }
    }else{
        matUR=this->udvRcrd.BUTZ[udvInd];
        matVR=this->udvRcrd.BVTZ[udvInd];
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            diagDR[siteInd]=this->udvRcrd.BDTZ[udvInd].coeffRef(siteInd,siteInd);
        }
    }
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        if(abs(diagDR[siteInd])>=1){
            diagDRLargeInv[siteInd]=1/diagDR[siteInd];
            diagDRSmall[siteInd]=1;
        }else{
            diagDRLargeInv[siteInd]=1;
            diagDRSmall[siteInd]=diagDR[siteInd];
        }
    }
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDRLargeInv[siteInd]);
    }
    matDRLargeInv.setFromTriplets(triList.begin(),triList.end());
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDRSmall[siteInd]);
    }
    matDRSmall.setFromTriplets(triList.begin(),triList.end());

    MatrixXd matUL(this->model.numSite,this->model.numSite);
    MatrixXd matVL(this->model.numSite,this->model.numSite);
    vector<double> diagDL(this->model.numSite);
    vector<double> diagDLLargeInv(this->model.numSite);
    vector<double> diagDLSmall(this->model.numSite);
    SparseMatrix<double> matDLLargeInv(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDLSmall(this->model.numSite,this->model.numSite);
    if(udvInd==this->numUdv-1){
        matUL.setIdentity();
        matVL.setIdentity();
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            diagDL[siteInd]=1;
        }
    }else{
        matUL=this->udvRcrd.BUBT[this->numUdv-2-udvInd];
        matVL=this->udvRcrd.BVBT[this->numUdv-2-udvInd];
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            diagDL[siteInd]=this->udvRcrd.BDBT[this->numUdv-2-udvInd].coeffRef(siteInd,siteInd);
        }
    }
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        if(abs(diagDL[siteInd])>=1){
            diagDLLargeInv[siteInd]=1/diagDL[siteInd];
            diagDLSmall[siteInd]=1;
        }else{
            diagDLLargeInv[siteInd]=1;
            diagDLSmall[siteInd]=diagDL[siteInd];
        }
    }
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDLLargeInv[siteInd]);
    }
    matDLLargeInv.setFromTriplets(triList.begin(),triList.end());
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDLSmall[siteInd]);
    }
    matDLSmall.setFromTriplets(triList.begin(),triList.end());
    
    MatrixXd matULInv=matUL.inverse();
    MatrixXd matURInv=matUR.inverse();
    MatrixXd matMid1,matMid2,matMid;
    matMid1=matURInv*matULInv;
    matMid1=matDRLargeInv*matMid1;
    matMid1=matMid1*matDLLargeInv;
    matMid2=matVR*matVL;
    matMid2=matDRSmall*matMid2;
    matMid2=matMid2*matDLSmall;
    matMid=matMid1+matMid2;
    matMid=matMid.inverse();
    this->grnFunc.equB=matMid;
    this->grnFunc.equB=matDLLargeInv*this->grnFunc.equB;
    this->grnFunc.equB=matULInv*this->grnFunc.equB;
    this->grnFunc.equB=this->grnFunc.equB*matDRLargeInv;
    this->grnFunc.equB=this->grnFunc.equB*matURInv;
}
void BscMC::syncGrnEqu(const int udvInd){
    this->syncGrnEquA(udvInd);
    this->syncGrnEquB(udvInd);
}

void BscMC::syncGrnModEquA(const int udvInd){
    this->syncGrnEquA(udvInd);
    this->grnFunc.modEquA=this->grnFunc.equA;
    this->grnFunc.modEquA=this->matK.A1PQuar*this->grnFunc.modEquA;
    this->grnFunc.modEquA=this->matK.A2PHalf*this->grnFunc.modEquA;
    this->grnFunc.modEquA=this->matK.A1PQuar*this->grnFunc.modEquA;
    this->grnFunc.modEquA=this->grnFunc.modEquA*this->matK.A1NQuar;
    this->grnFunc.modEquA=this->grnFunc.modEquA*this->matK.A2NHalf;
    this->grnFunc.modEquA=this->grnFunc.modEquA*this->matK.A1NQuar;
}
void BscMC::syncGrnModEquB(const int udvInd){
    this->syncGrnEquB(udvInd);
    this->grnFunc.modEquB=this->grnFunc.equB;
    this->grnFunc.modEquB=this->matK.B1PQuar*this->grnFunc.modEquB;
    this->grnFunc.modEquB=this->matK.B2PHalf*this->grnFunc.modEquB;
    this->grnFunc.modEquB=this->matK.B1PQuar*this->grnFunc.modEquB;
    this->grnFunc.modEquB=this->grnFunc.modEquB*this->matK.B1NQuar;
    this->grnFunc.modEquB=this->grnFunc.modEquB*this->matK.B2NHalf;
    this->grnFunc.modEquB=this->grnFunc.modEquB*this->matK.B1NQuar;
}
void BscMC::syncGrnModEqu(const int udvInd){
    this->syncGrnModEquA(udvInd);
    this->syncGrnModEquB(udvInd);
}


void BscMC::prepScanFromZero(){
    this->syncUdvRcrd();
    this->syncGrnModEqu(-1);
}
void BscMC::prepScanFromBeta(){
    this->syncUdvRcrd();
    this->syncGrnModEqu(this->numUdv-1);
}

/*calculate weight ratio of new and old configuration*/
/*this demands updated modified Green's functions*/
double BscMC::wgtRatFer(const int siteUpd,const double stateChng)const{
    double expTmpA=exp(-this->timeIntv*this->model.strCpl*stateChng);
    double expTmpB=1/expTmpA;
    double resA=(expTmpA-1)*(1-this->grnFunc.modEquA(siteUpd,siteUpd))+1;
    double resB=(expTmpB-1)*(1-this->grnFunc.modEquB(siteUpd,siteUpd))+1;
    double res=resA*resA*resB*resB;
    return res;
}
double BscMC::wgtRatBos(const int timeUpd,const int siteUpd,const double stateChng)const{
    const double stateOld=this->conf[timeUpd][siteUpd];
    const double stateNew=stateOld+stateChng;

    /*dynamic term*/
    double diffDyna=0;
    if(timeUpd==0){
        diffDyna=pow(stateNew,2)-pow(stateOld,2)-2*stateChng*conf[timeUpd+1][siteUpd];
    }
    if(timeUpd==(this->numSlc-1)){
        diffDyna=pow(stateNew,2)-pow(stateOld,2)-2*stateChng*conf[timeUpd-1][siteUpd];
    }
    if((timeUpd>0)&&(timeUpd<(this->numSlc-1))){
        diffDyna=2*pow(stateNew,2)-2*pow(stateOld,2);
        diffDyna=diffDyna-2*stateChng*conf[timeUpd-1][siteUpd];
        diffDyna=diffDyna-2*stateChng*conf[timeUpd+1][siteUpd];
    }

    /*gradient term*/
    int neigSite[4];
    this->model.genNeigSite(siteUpd,neigSite);
    double sumTmp=0;
    for(int i=0;i<4;i++){
        sumTmp+=conf[timeUpd][neigSite[i]];
    }
    double diffGrad=4*pow(stateNew,2)-4*pow(stateOld,2);
    diffGrad=diffGrad-2*stateChng*sumTmp;

    /*square and quartic term*/
    double diffSq=pow(stateNew,2)-pow(stateOld,2);
    double diffQu=0.5*pow(stateNew,4)-0.5*pow(stateOld,4);

    double tmp=0;
    tmp+=this->model.strDyn*diffDyna/this->timeIntv/this->timeIntv;
    tmp+=this->model.strGrad*diffGrad;
    tmp+=this->model.strSq*diffSq;
    tmp+=this->model.strQt*diffQu;
    return exp(-0.5*this->timeIntv*tmp);
}
double BscMC::wgtRat(const int timeUpd,const int siteUpd,const double stateChng)const{
    double resFer=this->wgtRatFer(siteUpd,stateChng);
    double resBos=this->wgtRatBos(timeUpd,siteUpd,stateChng);
    return resFer*resBos;
}

/*update modEquGre when a Metropolis proposal is accepted*/
void BscMC::updGrnModEqu(const int siteUpd,const double stateChng){
    double vChng,val;
    vector<double> matMulEffRow(this->model.numSite);
    vector< Triplet<double> > triList(this->model.numSite);
    SparseMatrix<double> matMul(this->model.numSite,this->model.numSite);

    /*for orbital A*/
    vChng=exp(-this->timeIntv*this->model.strCpl*stateChng)-1;
    val=vChng*(1-this->grnFunc.modEquA(siteUpd,siteUpd))+1;
    val=vChng/val;
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        matMulEffRow[siteInd]=this->grnFunc.modEquA(siteUpd,siteInd);
    }
    matMulEffRow[siteUpd]-=1;
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        triList[siteInd]=Triplet<double>(siteUpd,siteInd,val*matMulEffRow[siteInd]);
    }
    matMul.setFromTriplets(triList.begin(),triList.end());
    this->grnFunc.modEquA=this->grnFunc.modEquA+this->grnFunc.modEquA*matMul;

    /*for orbital B*/
    vChng=exp(this->timeIntv*this->model.strCpl*stateChng)-1;
    val=vChng*(1-this->grnFunc.modEquB(siteUpd,siteUpd))+1;
    val=vChng/val;
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        matMulEffRow[siteInd]=this->grnFunc.modEquB(siteUpd,siteInd);
    }
    matMulEffRow[siteUpd]-=1;
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        triList[siteInd]=Triplet<double>(siteUpd,siteInd,val*matMulEffRow[siteInd]);
    }
    matMul.setFromTriplets(triList.begin(),triList.end());
    this->grnFunc.modEquB=this->grnFunc.modEquB+this->grnFunc.modEquB*matMul;
}

void BscMC::swpSlc(const int slcUpd){
    double stateChng,accpProb;
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        stateChng=this->win*(2*this->ranMT.ranReal()-1);
        accpProb=this->wgtRat(slcUpd,siteInd,stateChng);
        if(this->ranMT.ranReal()<accpProb){
            this->conf[slcUpd][siteInd]+=stateChng;
            this->updGrnModEqu(siteInd,stateChng);
        }
    }
}
void BscMC::swpSlc(const int slcUpd,double& accpRate){
    double stateChng,accProb;
    int numAccp=0;
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        stateChng=win*(2*this->ranMT.ranReal()-1);
        accProb=this->wgtRat(slcUpd,siteInd,stateChng);
        if(this->ranMT.ranReal()<accProb){
            this->conf[slcUpd][siteInd]+=stateChng;
            this->updGrnModEqu(siteInd,stateChng);
            numAccp++;
        }
    }
    accpRate=double(numAccp)/this->model.numSite;
}

/*evolve Green's functions forward (slice index ++) or backward (slice index --)*/
void BscMC::evoGrnFwd(const int slcInd){
    this->grnFunc.modEquA=this->matK.A1NHalf*this->grnFunc.modEquA;
    this->grnFunc.modEquA=this->matK.A2N*this->grnFunc.modEquA;
    this->grnFunc.modEquA=this->matK.A1NHalf*this->grnFunc.modEquA;
    this->grnFunc.modEquA=this->genMatExpV(slcInd,0)*this->grnFunc.modEquA;

    this->grnFunc.modEquA=this->grnFunc.modEquA*this->matK.A1PHalf;
    this->grnFunc.modEquA=this->grnFunc.modEquA*this->matK.A2P;
    this->grnFunc.modEquA=this->grnFunc.modEquA*this->matK.A1PHalf;
    this->grnFunc.modEquA=this->grnFunc.modEquA*this->genMatExpV(slcInd,1);

    this->grnFunc.modEquB=this->matK.B1NHalf*this->grnFunc.modEquB;
    this->grnFunc.modEquB=this->matK.B2N*this->grnFunc.modEquB;
    this->grnFunc.modEquB=this->matK.B1NHalf*this->grnFunc.modEquB;
    this->grnFunc.modEquB=this->genMatExpV(slcInd,1)*this->grnFunc.modEquB;

    this->grnFunc.modEquB=this->grnFunc.modEquB*this->matK.B1PHalf;
    this->grnFunc.modEquB=this->grnFunc.modEquB*this->matK.B2P;
    this->grnFunc.modEquB=this->grnFunc.modEquB*this->matK.B1PHalf;
    this->grnFunc.modEquB=this->grnFunc.modEquB*this->genMatExpV(slcInd,0);
}
void BscMC::evoGrnBwd(const int slcInd){
    this->grnFunc.modEquA=this->genMatExpV(slcInd,1)*this->grnFunc.modEquA;
    this->grnFunc.modEquA=this->matK.A1PHalf*this->grnFunc.modEquA;
    this->grnFunc.modEquA=this->matK.A2P*this->grnFunc.modEquA;
    this->grnFunc.modEquA=this->matK.A1PHalf*this->grnFunc.modEquA;

    this->grnFunc.modEquA=this->grnFunc.modEquA*this->genMatExpV(slcInd,0);
    this->grnFunc.modEquA=this->grnFunc.modEquA*this->matK.A1NHalf;
    this->grnFunc.modEquA=this->grnFunc.modEquA*this->matK.A2N;
    this->grnFunc.modEquA=this->grnFunc.modEquA*this->matK.A1NHalf;

    this->grnFunc.modEquB=this->genMatExpV(slcInd,0)*this->grnFunc.modEquB;
    this->grnFunc.modEquB=this->matK.B1PHalf*this->grnFunc.modEquB;
    this->grnFunc.modEquB=this->matK.B2P*this->grnFunc.modEquB;
    this->grnFunc.modEquB=this->matK.B1PHalf*this->grnFunc.modEquB;

    this->grnFunc.modEquB=this->grnFunc.modEquB*this->genMatExpV(slcInd,1);
    this->grnFunc.modEquB=this->grnFunc.modEquB*this->matK.B1NHalf;
    this->grnFunc.modEquB=this->grnFunc.modEquB*this->matK.B2N;
    this->grnFunc.modEquB=this->grnFunc.modEquB*this->matK.B1NHalf;
}

void BscMC::scanZeroBeta(){
    MatrixXd prodA(this->model.numSite,this->model.numSite);
    MatrixXd prodB(this->model.numSite,this->model.numSite);
    prodA.setIdentity();
    prodB.setIdentity();

    MatrixXd matUA(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDA(this->model.numSite,this->model.numSite);
    MatrixXd matVA(this->model.numSite,this->model.numSite);
    matUA.setIdentity();
    matDA.setIdentity();
    matVA.setIdentity();

    MatrixXd matUB(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDB(this->model.numSite,this->model.numSite);
    MatrixXd matVB(this->model.numSite,this->model.numSite);
    matUB.setIdentity();
    matDB.setIdentity();
    matVB.setIdentity();

    int udvInd;
    MatrixXd matTmp(this->model.numSite,this->model.numSite);
    MatrixXd matUTmp(this->model.numSite,this->model.numSite);
    MatrixXd matVTmp(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDTmp(this->model.numSite,this->model.numSite);

    for(int slcInd=0;slcInd<this->numSlc;slcInd++){
        this->evoGrnFwd(slcInd);
        this->swpSlc(slcInd);

        
        prodA=matK.A1NQuar*prodA;
        prodA=matK.A2NHalf*prodA;
        prodA=matK.A1NQuar*prodA;
        prodA=genMatExpV(slcInd,0)*prodA;
        prodA=matK.A1NQuar*prodA;
        prodA=matK.A2NHalf*prodA;
        prodA=matK.A1NQuar*prodA;

        prodB=matK.B1NQuar*prodB;
        prodB=matK.B2NHalf*prodB;
        prodB=matK.B1NQuar*prodB;
        prodB=genMatExpV(slcInd,1)*prodB;
        prodB=matK.B1NQuar*prodB;
        prodB=matK.B2NHalf*prodB;
        prodB=matK.B1NQuar*prodB;
        
        if((slcInd+1)%this->udvIntv==0){
            udvInd=(slcInd+1)/this->udvIntv-1;

            /*update B-Matrice-Tau-Zero lists for orbital A*/
            matTmp=prodA*matUA;
            matTmp=matTmp*matDA;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matUA=matUTmp;
            matDA=matDTmp;
            matVA=matVTmp*matVA;
            this->udvRcrd.AUTZ[udvInd]=matUA;
            this->udvRcrd.ADTZ[udvInd]=matDA;
            this->udvRcrd.AVTZ[udvInd]=matVA;
            prodA.setIdentity();

            /*update B-Matrice-Tau-Zero lists for orbital B*/
            matTmp=prodB*matUB;
            matTmp=matTmp*matDB;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matUB=matUTmp;
            matDB=matDTmp;
            matVB=matVTmp*matVB;
            this->udvRcrd.BUTZ[udvInd]=matUB;
            this->udvRcrd.BDTZ[udvInd]=matDB;
            this->udvRcrd.BVTZ[udvInd]=matVB;
            prodB.setIdentity();

            /*recalculate Green's functions for orbital A and B*/
            this->syncGrnModEqu(udvInd);
        }
    }
}

void BscMC::scanBetaZero(){
    MatrixXd prodA(this->model.numSite,this->model.numSite);
    MatrixXd prodB(this->model.numSite,this->model.numSite);
    prodA.setIdentity();
    prodB.setIdentity();

    MatrixXd matUA(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDA(this->model.numSite,this->model.numSite);
    MatrixXd matVA(this->model.numSite,this->model.numSite);
    matUA.setIdentity();
    matDA.setIdentity();
    matVA.setIdentity();

    MatrixXd matUB(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDB(this->model.numSite,this->model.numSite);
    MatrixXd matVB(this->model.numSite,this->model.numSite);
    matUB.setIdentity();
    matDB.setIdentity();
    matVB.setIdentity();

    int udvInd;
    MatrixXd matTmp(this->model.numSite,this->model.numSite);
    MatrixXd matUTmp(this->model.numSite,this->model.numSite);
    MatrixXd matVTmp(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDTmp(this->model.numSite,this->model.numSite);

    for(int slcInd=this->numSlc-1;slcInd>=0;slcInd--){
        this->swpSlc(slcInd);

        prodA=matK.A1NQuar*prodA;
        prodA=matK.A2NHalf*prodA;
        prodA=matK.A1NQuar*prodA;
        prodA=genMatExpV(slcInd,0)*prodA;
        prodA=matK.A1NQuar*prodA;
        prodA=matK.A2NHalf*prodA;
        prodA=matK.A1NQuar*prodA;

        prodB=matK.B1NQuar*prodB;
        prodB=matK.B2NHalf*prodB;
        prodB=matK.B1NQuar*prodB;
        prodB=genMatExpV(slcInd,1)*prodB;
        prodB=matK.B1NQuar*prodB;
        prodB=matK.B2NHalf*prodB;
        prodB=matK.B1NQuar*prodB;

        if((this->numSlc-slcInd)%this->udvIntv==0){
            udvInd=(this->numSlc-slcInd)/this->udvIntv-1;

            matTmp=prodA*matUA;
            matTmp=matTmp*matDA;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matUA=matUTmp;
            matDA=matDTmp;
            matVA=matVTmp*matVA;
            this->udvRcrd.AUBT[udvInd]=matUA.transpose();
            this->udvRcrd.ADBT[udvInd]=matDA;
            this->udvRcrd.AVBT[udvInd]=matVA.transpose();
            prodA.setIdentity();

            matTmp=prodB*matUB;
            matTmp=matTmp*matDB;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matUB=matUTmp;
            matDB=matDTmp;
            matVB=matVTmp*matVB;
            this->udvRcrd.BUBT[udvInd]=matUB.transpose();
            this->udvRcrd.BDBT[udvInd]=matDB;
            this->udvRcrd.BVBT[udvInd]=matVB.transpose();
            prodB.setIdentity();
        }  

        if((this->numSlc-slcInd)%this->udvIntv==0){
            this->syncGrnModEqu(this->numUdv-2-udvInd);
        }else{
            this->evoGrnBwd(slcInd);
        }
    }
}

void BscMC::scanZeroBeta(double& accpRate,double& errGrn,double& condNum){
    double accpRateSum=0;
    double accpRateTmp;
    double errGrnSum=0;
    double condNumSum=0;


    MatrixXd prodA(this->model.numSite,this->model.numSite);
    MatrixXd prodB(this->model.numSite,this->model.numSite);
    prodA.setIdentity();
    prodB.setIdentity();

    MatrixXd matUA(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDA(this->model.numSite,this->model.numSite);
    MatrixXd matVA(this->model.numSite,this->model.numSite);
    matUA.setIdentity();
    matDA.setIdentity();
    matVA.setIdentity();

    MatrixXd matUB(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDB(this->model.numSite,this->model.numSite);
    MatrixXd matVB(this->model.numSite,this->model.numSite);
    matUB.setIdentity();
    matDB.setIdentity();
    matVB.setIdentity();

    int udvInd;
    MatrixXd matTmp(this->model.numSite,this->model.numSite);
    MatrixXd matUTmp(this->model.numSite,this->model.numSite);
    MatrixXd matVTmp(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDTmp(this->model.numSite,this->model.numSite);

    for(int slcInd=0;slcInd<this->numSlc;slcInd++){
        this->evoGrnFwd(slcInd);
        this->swpSlc(slcInd,accpRateTmp);
        accpRateSum+=accpRateTmp;
        
        prodA=matK.A1NQuar*prodA;
        prodA=matK.A2NHalf*prodA;
        prodA=matK.A1NQuar*prodA;
        prodA=this->genMatExpV(slcInd,0)*prodA;
        prodA=matK.A1NQuar*prodA;
        prodA=matK.A2NHalf*prodA;
        prodA=matK.A1NQuar*prodA;

        prodB=matK.B1NQuar*prodB;
        prodB=matK.B2NHalf*prodB;
        prodB=matK.B1NQuar*prodB;
        prodB=this->genMatExpV(slcInd,1)*prodB;
        prodB=matK.B1NQuar*prodB;
        prodB=matK.B2NHalf*prodB;
        prodB=matK.B1NQuar*prodB;
        
        if((slcInd+1)%this->udvIntv==0){
            /*accumulate effective conditional number*/
            JacobiSVD< MatrixXd > JacSVD;
            JacSVD.compute(prodA);
            VectorXd sinValA=JacSVD.singularValues();
            condNumSum+=sinValA(0)/sinValA(this->model.numSite-1);
            JacSVD.compute(prodB);
            VectorXd sinValB=JacSVD.singularValues();
            condNumSum+=sinValB(0)/sinValB(this->model.numSite-1);


            udvInd=(slcInd+1)/this->udvIntv-1;

            /*update B-Matrice-Tau-Zero lists for orbital A*/
            matTmp=prodA*matUA;
            matTmp=matTmp*matDA;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matUA=matUTmp;
            matDA=matDTmp;
            matVA=matVTmp*matVA;
            this->udvRcrd.AUTZ[udvInd]=matUA;
            this->udvRcrd.ADTZ[udvInd]=matDA;
            this->udvRcrd.AVTZ[udvInd]=matVA;
            prodA.setIdentity();

            /*update B-Matrice-Tau-Zero lists for orbital B*/
            matTmp=prodB*matUB;
            matTmp=matTmp*matDB;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matUB=matUTmp;
            matDB=matDTmp;
            matVB=matVTmp*matVB;
            this->udvRcrd.BUTZ[udvInd]=matUB;
            this->udvRcrd.BDTZ[udvInd]=matDB;
            this->udvRcrd.BVTZ[udvInd]=matVB;
            prodB.setIdentity();

            /*record error of modified Green's function*/
            MatrixXd diffModEquGreA=this->grnFunc.modEquA;
            MatrixXd diffModEquGreB=this->grnFunc.modEquB;

            /*recalculate Green's functions for orbital A and B*/
            this->syncGrnModEqu(udvInd);
            diffModEquGreA=diffModEquGreA-this->grnFunc.modEquA;
            diffModEquGreB=diffModEquGreB-this->grnFunc.modEquB;
            // double diffSumA=0;
            // double diffSumB=0;
            // for(int i=0;i<this->model.numSite;i++){
            //     for(int j=0;j<this->model.numSite;j++){
            //         diffSumA+=abs(diffModEquGreA(i,j));
            //         diffSumB+=abs(diffModEquGreB(i,j));
            //     }
            // }
            // errGrnSum+=diffSumA;
            // errGrnSum+=diffSumB;
            errGrnSum+=diffModEquGreA.norm();
            errGrnSum+=diffModEquGreB.norm();
        }
    }

    accpRate=accpRateSum/this->numSlc;
    condNum=condNumSum/2/this->numUdv;
    errGrn=errGrnSum/2/this->numUdv;
}

void BscMC::scanBetaZero(double& accpRate,double& errGrn,double& condNum){
    double accpRateSum=0;
    double accpRateTmp;
    double errGreSum=0;
    double condNumSum=0;

    MatrixXd prodA(this->model.numSite,this->model.numSite);
    MatrixXd prodB(this->model.numSite,this->model.numSite);
    prodA.setIdentity();
    prodB.setIdentity();

    MatrixXd matUA(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDA(this->model.numSite,this->model.numSite);
    MatrixXd matVA(this->model.numSite,this->model.numSite);
    matUA.setIdentity();
    matDA.setIdentity();
    matVA.setIdentity();

    MatrixXd matUB(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDB(this->model.numSite,this->model.numSite);
    MatrixXd matVB(this->model.numSite,this->model.numSite);
    matUB.setIdentity();
    matDB.setIdentity();
    matVB.setIdentity();

    int udvInd;
    MatrixXd matTmp(this->model.numSite,this->model.numSite);
    MatrixXd matUTmp(this->model.numSite,this->model.numSite);
    MatrixXd matVTmp(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDTmp(this->model.numSite,this->model.numSite);

    for(int slcInd=this->numSlc-1;slcInd>=0;slcInd--){
        this->swpSlc(slcInd,accpRateTmp);
        accpRateSum+=accpRateTmp;

        prodA=matK.A1NQuar*prodA;
        prodA=matK.A2NHalf*prodA;
        prodA=matK.A1NQuar*prodA;
        prodA=this->genMatExpV(slcInd,0)*prodA;
        prodA=matK.A1NQuar*prodA;
        prodA=matK.A2NHalf*prodA;
        prodA=matK.A1NQuar*prodA;

        prodB=matK.B1NQuar*prodB;
        prodB=matK.B2NHalf*prodB;
        prodB=matK.B1NQuar*prodB;
        prodB=this->genMatExpV(slcInd,1)*prodB;
        prodB=matK.B1NQuar*prodB;
        prodB=matK.B2NHalf*prodB;
        prodB=matK.B1NQuar*prodB;

        if((this->numSlc-slcInd)%this->udvIntv==0){
            /*accumulate effective conditional number*/
            JacobiSVD< MatrixXd > JacSVD;
            JacSVD.compute(prodA);
            VectorXd sinValA=JacSVD.singularValues();
            condNumSum+=sinValA(0)/sinValA(this->model.numSite-1);
            JacSVD.compute(prodB);
            VectorXd sinValB=JacSVD.singularValues();
            condNumSum+=sinValB(0)/sinValB(this->model.numSite-1);

            udvInd=(this->numSlc-slcInd)/this->udvIntv-1;

            matTmp=prodA*matUA;
            matTmp=matTmp*matDA;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matUA=matUTmp;
            matDA=matDTmp;
            matVA=matVTmp*matVA;
            this->udvRcrd.AUBT[udvInd]=matUA.transpose();
            this->udvRcrd.ADBT[udvInd]=matDA;
            this->udvRcrd.AVBT[udvInd]=matVA.transpose();
            prodA.setIdentity();

            matTmp=prodB*matUB;
            matTmp=matTmp*matDB;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matUB=matUTmp;
            matDB=matDTmp;
            matVB=matVTmp*matVB;
            this->udvRcrd.BUBT[udvInd]=matUB.transpose();
            this->udvRcrd.BDBT[udvInd]=matDB;
            this->udvRcrd.BVBT[udvInd]=matVB.transpose();
            prodB.setIdentity();
        }  

        if((this->numSlc-slcInd)%this->udvIntv==0){
            /*record error of modified Green's function*/
            this->evoGrnBwd(slcInd);
            MatrixXd diffModEquGreA=this->grnFunc.modEquA;
            MatrixXd diffModEquGreB=this->grnFunc.modEquB;
            
            this->syncGrnModEqu(this->numUdv-2-udvInd);

            diffModEquGreA=diffModEquGreA-this->grnFunc.modEquA;
            diffModEquGreB=diffModEquGreB-this->grnFunc.modEquB;
            // double diffSumA=0;
            // double diffSumB=0;
            // for(int i=0;i<this->model.numSite;i++){
            //     for(int j=0;j<this->model.numSite;j++){
            //         diffSumA+=abs(diffModEquGreA(i,j));
            //         diffSumB+=abs(diffModEquGreB(i,j));
            //     }
            // }
            // errGreSum+=diffSumA;
            // errGreSum+=diffSumB;
            errGreSum+=diffModEquGreA.norm();
            errGreSum+=diffModEquGreB.norm();
        }else{
            this->evoGrnBwd(slcInd);
        }
    }
    accpRate=accpRateSum/this->numSlc;
    errGrn=errGreSum/2/this->numUdv;
    condNum=condNumSum/2/this->numUdv;
}

void BscMC::cycScan(const int dir){
    if(dir==0){
        this->scanZeroBeta();
        this->scanBetaZero();
    }else{
        this->scanBetaZero();
        this->scanZeroBeta();
    }
}
void BscMC::cycScan(double& accpRate,double& errGrn,double& condNum,const int dir){
    double accpRateTmp;
    double errGrnTmp;
    double condNumTmp;
    if(dir==0){
        this->scanZeroBeta(accpRateTmp,errGrnTmp,condNumTmp);
        accpRate=accpRateTmp;
        errGrn=errGrnTmp;
        condNum=condNumTmp;

        this->scanBetaZero(accpRateTmp,errGrnTmp,condNumTmp);
        accpRate+=accpRateTmp;
        errGrn+=errGrnTmp;
        condNum+=condNumTmp;
    }else{
        this->scanBetaZero(accpRateTmp,errGrnTmp,condNumTmp);
        accpRate=accpRateTmp;
        errGrn=errGrnTmp;
        condNum=condNumTmp;

        this->scanZeroBeta(accpRateTmp,errGrnTmp,condNumTmp);
        accpRate+=accpRateTmp;
        errGrn+=errGrnTmp;
        condNum+=condNumTmp;
    }
    accpRate=accpRate/2;
    errGrn=errGrn/2;
    condNum=condNum/2;
}

void BscMC::printRoundOffInfo(ostream& iout){
    this->prepScanFromZero();
    double accpRate,errGrn,condNum;
    this->cycScan(accpRate,errGrn,condNum,0);
    iout<<"accumulate error of Green's functions: "<<errGrn<<endl;
    iout<<"effective UDV conditional number:      "<<condNum<<endl;
}

double BscMC::calcAccpRate(){
    this->prepScanFromZero();
    vector<double> arPath(0);
    double arTmp,bTmp,cTmp;
    for(int genInd=0;genInd<100;genInd++){
        this->cycScan(arTmp,bTmp,cTmp,0);
        arPath.push_back(arTmp);
        if((genInd+1)>=20&&(genInd+1)%5==0&&calcErrBin(arPath,5)<0.01){
            break;
        }
    }
    return calcAvrg(arPath);
}

void BscMC::adjWin(){
    double accpRate,bTmp,cTmp;
    this->cycScan(accpRate,bTmp,cTmp);
    if(accpRate<0.5){
        this->win=this->win/1.05;
    }else{
        this->win=this->win*1.05;
    }
}
void BscMC::detWin(const double winMin,const double winMax,ostream& iout){
    double winL=winMin;
    double winR=winMax;
    double winMid;
    double accpRateCurr;
    
    win=winL;
    const double accRateL=this->calcAccpRate();
    win=winR;
    const double accRateR=this->calcAccpRate();
    iout<<"accept rate for the left bound : "<<accRateL<<endl;
    iout<<"accept rate for the right bound: "<<accRateR<<endl;
    if(accRateL<0.45){
        iout<<"NOT FOUND: decrease the left bound!"<<endl;
        return;
    }
    if(accRateR>0.49){
        iout<<"NOT FOUND: increase the right bound!"<<endl;
        return;
    }

    for(int genInd=0;genInd<1000;genInd++){
        winMid=(winL+winR)/2;
        win=winMid;
        accpRateCurr=this->calcAccpRate();
        if(accpRateCurr<=0.45){
            winR=winMid;
        }
        if(accpRateCurr>=0.49){
            winL=winMid;
        }
        iout<<"trial window: "<<win<<"     accept rate: "<<accpRateCurr<<endl;
        if(accpRateCurr>0.45&&accpRateCurr<0.49){
            break;
        }
    }
    iout<<"Random window is fixed at "<<win<<" with average accept rate "<<accpRateCurr<<"."<<endl;

}

void BscMC::mesrEquGrn(MatrixXd* resA,MatrixXd* resB){
    MatrixXd matTmp(this->model.numSite,this->model.numSite);
    MatrixXd matUTmp(this->model.numSite,this->model.numSite);
    MatrixXd matVTmp(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDTmp(this->model.numSite,this->model.numSite);

    MatrixXd prodA(this->model.numSite,this->model.numSite);
    MatrixXd prodB(this->model.numSite,this->model.numSite);

    MatrixXd matUA(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDA(this->model.numSite,this->model.numSite);
    MatrixXd matVA(this->model.numSite,this->model.numSite);
    MatrixXd matUB(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDB(this->model.numSite,this->model.numSite);
    MatrixXd matVB(this->model.numSite,this->model.numSite);

    /*******************calculate B(\tau_l,0) for all slice*******************/
    prodA.setIdentity();
    prodB.setIdentity();

    matUA.setIdentity();
    matDA.setIdentity();
    matVA.setIdentity();
    
    matUB.setIdentity();
    matDB.setIdentity();
    matVB.setIdentity();

    MatrixXd* matFwdUA=new MatrixXd[this->numSlc];
    MatrixXd* matFwdVA=new MatrixXd[this->numSlc];
    MatrixXd* matFwdUB=new MatrixXd[this->numSlc];
    MatrixXd* matFwdVB=new MatrixXd[this->numSlc];
    SparseMatrix<double>* matFwdDA=new SparseMatrix<double>[this->numSlc];
    SparseMatrix<double>* matFwdDB=new SparseMatrix<double>[this->numSlc];

    for(int slcInd=0;slcInd<this->numSlc;slcInd++){
        prodA=matK.A1NQuar*prodA;
        prodA=matK.A2NHalf*prodA;
        prodA=matK.A1NQuar*prodA;
        prodA=genMatExpV(slcInd,0)*prodA;
        prodA=matK.A1NQuar*prodA;
        prodA=matK.A2NHalf*prodA;
        prodA=matK.A1NQuar*prodA;

        prodB=matK.B1NQuar*prodB;
        prodB=matK.B2NHalf*prodB;
        prodB=matK.B1NQuar*prodB;
        prodB=genMatExpV(slcInd,1)*prodB;
        prodB=matK.B1NQuar*prodB;
        prodB=matK.B2NHalf*prodB;
        prodB=matK.B1NQuar*prodB;

        matFwdUA[slcInd]=prodA*matUA;
        matFwdDA[slcInd]=matDA;
        matFwdVA[slcInd]=matVA;
        matFwdUB[slcInd]=prodB*matUB;
        matFwdDB[slcInd]=matDB;
        matFwdVB[slcInd]=matVB;
        
        if((slcInd+1)%this->udvIntv==0){
            /*update B-Matrice-Tau-Zero lists for orbital A*/
            matTmp=prodA*matUA;
            matTmp=matTmp*matDA;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matUA=matUTmp;
            matDA=matDTmp;
            matVA=matVTmp*matVA;
            prodA.setIdentity();

            /*update B-Matrice-Tau-Zero lists for orbital B*/
            matTmp=prodB*matUB;
            matTmp=matTmp*matDB;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matUB=matUTmp;
            matDB=matDTmp;
            matVB=matVTmp*matVB;
            prodB.setIdentity();
        }
    }

    /*******************calculate B(\beta,\tau_l) for all slice*******************/
    prodA.setIdentity();
    prodB.setIdentity();

    matUA.setIdentity();
    matDA.setIdentity();
    matVA.setIdentity();

    matUB.setIdentity();
    matDB.setIdentity();
    matVB.setIdentity();

    MatrixXd* matBwdUA=new MatrixXd[this->numSlc];
    MatrixXd* matBwdVA=new MatrixXd[this->numSlc];
    MatrixXd* matBwdUB=new MatrixXd[this->numSlc];
    MatrixXd* matBwdVB=new MatrixXd[this->numSlc];
    SparseMatrix<double>* matBwdDA=new SparseMatrix<double>[this->numSlc];
    SparseMatrix<double>* matBwdDB=new SparseMatrix<double>[this->numSlc];
    
    for(int slcInd=this->numSlc-1;slcInd>=0;slcInd--){
        prodA=matK.A1NQuar*prodA;
        prodA=matK.A2NHalf*prodA;
        prodA=matK.A1NQuar*prodA;
        prodA=genMatExpV(slcInd,0)*prodA;
        prodA=matK.A1NQuar*prodA;
        prodA=matK.A2NHalf*prodA;
        prodA=matK.A1NQuar*prodA;

        prodB=matK.B1NQuar*prodB;
        prodB=matK.B2NHalf*prodB;
        prodB=matK.B1NQuar*prodB;
        prodB=genMatExpV(slcInd,1)*prodB;
        prodB=matK.B1NQuar*prodB;
        prodB=matK.B2NHalf*prodB;
        prodB=matK.B1NQuar*prodB;

        matBwdUA[this->numSlc-1-slcInd]=(prodA*matUA).transpose();
        matBwdDA[this->numSlc-1-slcInd]=matDA;
        matBwdVA[this->numSlc-1-slcInd]=matVA.transpose();

        matBwdUB[this->numSlc-1-slcInd]=(prodB*matUB).transpose();
        matBwdDB[this->numSlc-1-slcInd]=matDB;
        matBwdVB[this->numSlc-1-slcInd]=matVB.transpose();

        if((this->numSlc-slcInd)%this->udvIntv==0){
            matTmp=prodA*matUA;
            matTmp=matTmp*matDA;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matUA=matUTmp;
            matDA=matDTmp;
            matVA=matVTmp*matVA;
            prodA.setIdentity();

            matTmp=prodB*matUB;
            matTmp=matTmp*matDB;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matUB=matUTmp;
            matDB=matDTmp;
            matVB=matVTmp*matVB;
            prodB.setIdentity();
        }  
    }

    /*==================calculate equal-time Green's function==============*/
    vector< Triplet<double> > triList(this->model.numSite);
    MatrixXd matUR(this->model.numSite,this->model.numSite);
    MatrixXd matVR(this->model.numSite,this->model.numSite);
    vector<double> diagDR(this->model.numSite);
    vector<double> diagDRLargeInv(this->model.numSite);
    vector<double> diagDRSmall(this->model.numSite);
    SparseMatrix<double> matDRLargeInv(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDRSmall(this->model.numSite,this->model.numSite);

    MatrixXd matUL(this->model.numSite,this->model.numSite);
    MatrixXd matVL(this->model.numSite,this->model.numSite);
    vector<double> diagDL(this->model.numSite);
    vector<double> diagDLLargeInv(this->model.numSite);
    vector<double> diagDLSmall(this->model.numSite);
    SparseMatrix<double> matDLLargeInv(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDLSmall(this->model.numSite,this->model.numSite);

    MatrixXd matULInv(this->model.numSite,this->model.numSite);
    MatrixXd matURInv(this->model.numSite,this->model.numSite);
    MatrixXd matMid(this->model.numSite,this->model.numSite);
    MatrixXd matMid1(this->model.numSite,this->model.numSite);
    MatrixXd matMid2(this->model.numSite,this->model.numSite);

    /*----------------for orbital A---------------*/
    for(int slcInd=0;slcInd<this->numSlc;slcInd++){
        matUR=matFwdUA[slcInd];
        matVR=matFwdVA[slcInd];
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            diagDR[siteInd]=matFwdDA[slcInd].coeffRef(siteInd,siteInd);
        }
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            if(abs(diagDR[siteInd])>=1){
                diagDRLargeInv[siteInd]=1/diagDR[siteInd];
                diagDRSmall[siteInd]=1;
            }else{
                diagDRLargeInv[siteInd]=1;
                diagDRSmall[siteInd]=diagDR[siteInd];
            }
        }
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDRLargeInv[siteInd]);
        }
        matDRLargeInv.setFromTriplets(triList.begin(),triList.end());
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDRSmall[siteInd]);
        }
        matDRSmall.setFromTriplets(triList.begin(),triList.end());

        if(slcInd==this->numSlc-1){
            matUL.setIdentity();
            matVL.setIdentity();
            for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
                diagDL[siteInd]=1;
            }
        }else{
            matUL=matBwdUA[this->numSlc-2-slcInd];
            matVL=matBwdVA[this->numSlc-2-slcInd];
            for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
                diagDL[siteInd]=matBwdDA[this->numSlc-2-slcInd].coeffRef(siteInd,siteInd);
            }
        }
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            if(abs(diagDL[siteInd])>=1){
                diagDLLargeInv[siteInd]=1/diagDL[siteInd];
                diagDLSmall[siteInd]=1;
            }else{
                diagDLLargeInv[siteInd]=1;
                diagDLSmall[siteInd]=diagDL[siteInd];
            }
        }
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDLLargeInv[siteInd]);
        }
        matDLLargeInv.setFromTriplets(triList.begin(),triList.end());
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDLSmall[siteInd]);
        }
        matDLSmall.setFromTriplets(triList.begin(),triList.end());

        matULInv=matUL.inverse();
        matURInv=matUR.inverse();
        matMid1=matURInv*matULInv;
        matMid1=matDRLargeInv*matMid1;
        matMid1=matMid1*matDLLargeInv;
        matMid2=matVR*matVL;
        matMid2=matDRSmall*matMid2;
        matMid2=matMid2*matDLSmall;
        matMid=matMid1+matMid2;
        matMid=matMid.inverse();
        resA[slcInd]=matMid;
        resA[slcInd]=matDLLargeInv*resA[slcInd];
        resA[slcInd]=matULInv*resA[slcInd];
        resA[slcInd]=resA[slcInd]*matDRLargeInv;
        resA[slcInd]=resA[slcInd]*matURInv;
    }

    /*----------------for orbital B---------------*/
    for(int slcInd=0;slcInd<this->numSlc;slcInd++){
        matUR=matFwdUB[slcInd];
        matVR=matFwdVB[slcInd];
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            diagDR[siteInd]=matFwdDB[slcInd].coeffRef(siteInd,siteInd);
        }
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            if(abs(diagDR[siteInd])>=1){
                diagDRLargeInv[siteInd]=1/diagDR[siteInd];
                diagDRSmall[siteInd]=1;
            }else{
                diagDRLargeInv[siteInd]=1;
                diagDRSmall[siteInd]=diagDR[siteInd];
            }
        }
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDRLargeInv[siteInd]);
        }
        matDRLargeInv.setFromTriplets(triList.begin(),triList.end());
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDRSmall[siteInd]);
        }
        matDRSmall.setFromTriplets(triList.begin(),triList.end());

        if(slcInd==this->numSlc-1){
            matUL.setIdentity();
            matVL.setIdentity();
            for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
                diagDL[siteInd]=1;
            }
        }else{
            matUL=matBwdUB[this->numSlc-2-slcInd];
            matVL=matBwdVB[this->numSlc-2-slcInd];
            for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
                diagDL[siteInd]=matBwdDB[this->numSlc-2-slcInd].coeffRef(siteInd,siteInd);
            }
        }
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            if(abs(diagDL[siteInd])>=1){
                diagDLLargeInv[siteInd]=1/diagDL[siteInd];
                diagDLSmall[siteInd]=1;
            }else{
                diagDLLargeInv[siteInd]=1;
                diagDLSmall[siteInd]=diagDL[siteInd];
            }
        }
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDLLargeInv[siteInd]);
        }
        matDLLargeInv.setFromTriplets(triList.begin(),triList.end());
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDLSmall[siteInd]);
        }
        matDLSmall.setFromTriplets(triList.begin(),triList.end());

        matULInv=matUL.inverse();
        matURInv=matUR.inverse();
        matMid1=matURInv*matULInv;
        matMid1=matDRLargeInv*matMid1;
        matMid1=matMid1*matDLLargeInv;
        matMid2=matVR*matVL;
        matMid2=matDRSmall*matMid2;
        matMid2=matMid2*matDLSmall;
        matMid=matMid1+matMid2;
        matMid=matMid.inverse();
        resB[slcInd]=matMid;
        resB[slcInd]=matDLLargeInv*resB[slcInd];
        resB[slcInd]=matULInv*resB[slcInd];
        resB[slcInd]=resB[slcInd]*matDRLargeInv;
        resB[slcInd]=resB[slcInd]*matURInv;
    }


    /*deallocate stuff newed in this sub-routine*/
    delete[] matFwdUA;
    delete[] matFwdDA;
    delete[] matFwdVA;
    delete[] matFwdUB;
    delete[] matFwdDB;
    delete[] matFwdVB;

    delete[] matBwdUA;
    delete[] matBwdDA;
    delete[] matBwdVA;
    delete[] matBwdUB;
    delete[] matBwdDB;
    delete[] matBwdVB;
}

void BscMC::mesrNeqGrn(MatrixXd* resA,MatrixXd* resB){
    /***************calculate B(\tau,0) for all slice****************/
    MatrixXd matTmp(this->model.numSite,this->model.numSite);
    MatrixXd matUTmp(this->model.numSite,this->model.numSite);
    MatrixXd matVTmp(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDTmp(this->model.numSite,this->model.numSite);

    MatrixXd prodA(this->model.numSite,this->model.numSite);
    MatrixXd prodB(this->model.numSite,this->model.numSite);

    MatrixXd matUA(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDA(this->model.numSite,this->model.numSite);
    MatrixXd matVA(this->model.numSite,this->model.numSite);
    MatrixXd matUB(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDB(this->model.numSite,this->model.numSite);
    MatrixXd matVB(this->model.numSite,this->model.numSite);

    prodA.setIdentity();
    prodB.setIdentity();

    matUA.setIdentity();
    matDA.setIdentity();
    matVA.setIdentity();
    
    matUB.setIdentity();
    matDB.setIdentity();
    matVB.setIdentity();

    MatrixXd* matFwdA=new MatrixXd[this->numSlc];
    MatrixXd* matFwdB=new MatrixXd[this->numSlc];

    for(int slcInd=0;slcInd<this->numSlc;slcInd++){
        prodA=matK.A1NQuar*prodA;
        prodA=matK.A2NHalf*prodA;
        prodA=matK.A1NQuar*prodA;
        prodA=genMatExpV(slcInd,0)*prodA;
        prodA=matK.A1NQuar*prodA;
        prodA=matK.A2NHalf*prodA;
        prodA=matK.A1NQuar*prodA;

        prodB=matK.B1NQuar*prodB;
        prodB=matK.B2NHalf*prodB;
        prodB=matK.B1NQuar*prodB;
        prodB=genMatExpV(slcInd,1)*prodB;
        prodB=matK.B1NQuar*prodB;
        prodB=matK.B2NHalf*prodB;
        prodB=matK.B1NQuar*prodB;

        matFwdA[slcInd]=prodA*matUA*matDA*matVA;
        matFwdB[slcInd]=prodB*matUB*matDB*matVB;
        
        if((slcInd+1)%this->udvIntv==0){
            /*update B-Matrice-Tau-Zero lists for orbital A*/
            matTmp=prodA*matUA;
            matTmp=matTmp*matDA;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matUA=matUTmp;
            matDA=matDTmp;
            matVA=matVTmp*matVA;
            prodA.setIdentity();

            /*update B-Matrice-Tau-Zero lists for orbital B*/
            matTmp=prodB*matUB;
            matTmp=matTmp*matDB;
            calcUdv(matTmp,matUTmp,matDTmp,matVTmp);
            matUB=matUTmp;
            matDB=matDTmp;
            matVB=matVTmp*matVB;
            prodB.setIdentity();
        }
    }

    /*calculate G(0) for orbital A and B*/
    MatrixXd equGrnZeroA(this->model.numSite,this->model.numSite);
    MatrixXd equGrnZeroB(this->model.numSite,this->model.numSite);
    vector<double> diagDA(this->model.numSite);
    vector<double> diagDB(this->model.numSite);
    vector<double> diagDALargeInv(this->model.numSite);
    vector<double> diagDBLargeInv(this->model.numSite);
    vector<double> diagDASmall(this->model.numSite);
    vector<double> diagDBSmall(this->model.numSite);
    SparseMatrix<double> matDALargeInv(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDBLargeInv(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDASmall(this->model.numSite,this->model.numSite);
    SparseMatrix<double> matDBSmall(this->model.numSite,this->model.numSite);

    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        diagDA[siteInd]=matDA.coeffRef(siteInd,siteInd);
        diagDB[siteInd]=matDB.coeffRef(siteInd,siteInd);
    }
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        if(abs(diagDA[siteInd])>=1){
            diagDALargeInv[siteInd]=1/diagDA[siteInd];
            diagDASmall[siteInd]=1;
        }else{
            diagDALargeInv[siteInd]=1;
            diagDASmall[siteInd]=diagDA[siteInd];
        }
        if(abs(diagDB[siteInd])>=1){
            diagDBLargeInv[siteInd]=1/diagDB[siteInd];
            diagDBSmall[siteInd]=1;
        }else{
            diagDBLargeInv[siteInd]=1;
            diagDBSmall[siteInd]=diagDB[siteInd];
        }
    }
    vector< Triplet<double> > triList(this->model.numSite);
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDALargeInv[siteInd]);
    }
    matDALargeInv.setFromTriplets(triList.begin(),triList.end());
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDASmall[siteInd]);
    }
    matDASmall.setFromTriplets(triList.begin(),triList.end());
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDBLargeInv[siteInd]);
    }
    matDBLargeInv.setFromTriplets(triList.begin(),triList.end());
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
        triList[siteInd]=Triplet<double>(siteInd,siteInd,diagDBSmall[siteInd]);
    }
    matDBSmall.setFromTriplets(triList.begin(),triList.end());

    MatrixXd matUAInv=matUA.inverse();
    MatrixXd matUBInv=matUB.inverse();
    equGrnZeroA=(matDALargeInv*matUAInv+matDASmall*matVA).inverse();
    equGrnZeroA=equGrnZeroA*matDALargeInv*matUAInv;
    equGrnZeroB=(matDBLargeInv*matUBInv+matDBSmall*matVB).inverse();
    equGrnZeroB=equGrnZeroB*matDBLargeInv*matUBInv;

    /*calculate time-displaced Green's function*/
    for(int slcInd=0;slcInd<this->numSlc-1;slcInd++){
        resA[slcInd+1]=matFwdA[slcInd]*equGrnZeroA;
        resB[slcInd+1]=matFwdB[slcInd]*equGrnZeroB;
    }
    resA[0]=equGrnZeroA;
    resB[0]=equGrnZeroB;

    /*deallocate stuff newed in this sub-routine*/
    delete[] matFwdA;
    delete[] matFwdB;
}

void BscMC::findSpacRSite(int rInd,int* res)const{
    int rY=rInd/this->model.lenX;
    int rX=rInd-rY*this->model.lenX;
    for(int iX=0;iX<this->model.lenX;iX++){
        for(int iY=0;iY<this->model.lenY;iY++){
            int jX=(rX+iX)%this->model.lenX;
            int jY=(rY+iY)%this->model.lenY;
            int i=iX+iY*this->model.lenX;
            int j=jX+jY*this->model.lenY;
            res[i]=j;
        }
    }
}

void BscMC::calcEquCorrSngSlc(MatrixXd& equGrnAIn,MatrixXd& equGrnBIn,vector<double>& res)const{
    res.resize(this->model.numSite);
    int siteSpacR[this->model.numSite];
    for(int rInd=0;rInd<this->model.numSite;rInd++){
        double tmp=0;
        if(rInd==0){
            for(int i=0;i<this->model.numSite;i++){
                tmp+=2*equGrnAIn(i,i)*(1+equGrnAIn(i,i));
                tmp+=2*equGrnBIn(i,i)*(1+equGrnBIn(i,i));
                tmp-=8*equGrnAIn(i,i)*equGrnBIn(i,i);
            }
            tmp=tmp/this->model.numSite;
        }else{
            this->findSpacRSite(rInd,siteSpacR);
            int j;
            for(int i=0;i<this->model.numSite;i++){
                j=siteSpacR[i];
                tmp+=4*equGrnAIn(i,i)*equGrnAIn(j,j);
                tmp-=2*equGrnAIn(i,j)*equGrnAIn(j,i);
                tmp+=4*equGrnBIn(i,i)*equGrnBIn(j,j);
                tmp-=2*equGrnBIn(i,j)*equGrnBIn(j,i);
                tmp-=4*equGrnAIn(i,i)*equGrnBIn(j,j);
                tmp-=4*equGrnAIn(j,j)*equGrnBIn(i,i);
            }
            tmp=tmp/this->model.numSite;
        }
        res[rInd]=tmp;
    }
}

void BscMC::calcEquCorr(MatrixXd* equGrnAIn,MatrixXd* equGrnBIn,vector<double>& res)const{
    res.resize(this->model.numSite);
    for(int i=0;i<res.size();i++){
        res[i]=0;
    }
    vector<double> resTmp;
    for(int slcInd=0;slcInd<this->numSlc;slcInd++){
        this->calcEquCorrSngSlc(equGrnAIn[slcInd],equGrnBIn[slcInd],resTmp);
        for(int i=0;i<this->model.numSite;i++){
            res[i]+=resTmp[i];
        }
    }

    for(int i=0;i<this->model.numSite;i++){
        res[i]=res[i]/this->numSlc;
    }
}

double BscMC::mesrConjFld()const{
    double res=0;
    for(int i=0;i<this->numSlc;i++){
        for(int j=0;j<this->model.numSite;j++){
            res+=this->conf[i][j];
        }
    }
    res=res/this->numSlc/this->model.numSite;
    return res;
}

double BscMC::mesrConjFldSq()const{
    double res=0;
    for(int i=0;i<this->numSlc;i++){
        double resSlc=0;
        for(int j=0;j<this->model.numSite;j++){
            resSlc+=this->conf[i][j];
        }
        resSlc=resSlc/this->model.numSite;
        res+=resSlc*resSlc;
    }
    return res/this->numSlc;
}

double BscMC::mesrConjFldQt()const{
    double res=0;
    for(int i=0;i<this->numSlc;i++){
        double resSlc=0;
        for(int j=0;j<this->model.numSite;j++){
            resSlc+=this->conf[i][j];
        }
        resSlc=resSlc/this->model.numSite;
        res+=pow(resSlc,4);
    }
    return res/this->numSlc;
}

double BscMC::calcEquCorrFMSngSlc(MatrixXd& equGrnAIn,MatrixXd& equGrnBIn)const{
    double trGA=equGrnAIn.trace();
    double trGB=equGrnBIn.trace();
    double trGASq=(equGrnAIn*equGrnAIn).trace();
    double trGBSq=(equGrnBIn*equGrnBIn).trace();
    double res=4*(trGA-trGB)*(trGA-trGB)+2*(trGA+trGB)-2*trGASq-2*trGBSq;
    return res/this->model.numSite/this->model.numSite;
}

double BscMC::calcEquCorrFM(MatrixXd* equGrnAIn,MatrixXd* equGrnBIn)const{
    vector<double> vecTmp(this->numSlc);
    for(int i=0;i<this->numSlc;i++){
        vecTmp[i]=this->calcEquCorrFMSngSlc(equGrnAIn[i],equGrnBIn[i]);
    }
    return calcAvrg(vecTmp);
}

double BscMC::calcOrdFMSngSlc(MatrixXd& equGrnAIn,MatrixXd& equGrnBIn)const{
    double res=equGrnBIn.trace()-equGrnAIn.trace();
    return res*2/this->model.numSite;
}
double BscMC::calcOrdFM(MatrixXd* equGrnAIn,MatrixXd* equGrnBIn)const{
    vector<double> vecTmp(this->numSlc);
    for(int i=0;i<this->numSlc;i++){
        vecTmp[i]=this->calcOrdFMSngSlc(equGrnAIn[i],equGrnBIn[i]);
    }
    return calcAvrg(vecTmp);
}
double BscMC::calcOrdFMSq(MatrixXd* equGrnAIn,MatrixXd* equGrnBIn)const{
    vector<double> vecTmp;
    double tmp;
    for(int i=0;i<this->numSlc;i++){
        tmp=this->calcOrdFMSngSlc(equGrnAIn[i],equGrnBIn[i]);
        vecTmp[i]=tmp*tmp;
    }
    return calcAvrg(vecTmp);
}
double BscMC::calcOrdFMQt(MatrixXd* equGrnAIn,MatrixXd* equGrnBIn)const{
    vector<double> vecTmp;
    double tmp;
    for(int i=0;i<this->numSlc;i++){
        tmp=this->calcOrdFMSngSlc(equGrnAIn[i],equGrnBIn[i]);
        vecTmp[i]=pow(tmp,4);
    }
    return calcAvrg(vecTmp);
}

void BscMC::exportStts(string loct){
    ofstream sout;
    sout.open(loct);
    sout<<this->model.lenX<<" "<<this->model.lenY<<endl;
    sout<<this->invTem<<endl;
    sout<<this->numUdv<<" "<<this->udvIntv<<endl;
    sout<<this->ranMT.readSeed()<<" "<<this->win<<endl;
    sout<<this->model.strHop<<" "<<this->model.ratHop<<endl;
    sout<<this->model.strCpl<<endl;
    sout<<this->model.strDyn<<" "<<this->model.strGrad<<" ";
    sout<<this->model.strSq<<" "<<this->model.strQt<<endl;
    sout<<this->model.chem<<" "<<this->model.biasChem<<endl;
    for(int i=0;i<this->numSlc;i++){
        for(int j=0;j<this->model.numSite;j++){
            sout<<this->conf[i][j]<<" ";
        }
        sout<<endl;
    }
    sout<<endl;
    sout.close();
}
void BscMC::importStts(string loct){
    ifstream sin;
    sin.open(loct);

    int pLenX,pLenY;
    sin>>pLenX>>pLenY;
    this->setLattSize(pLenX,pLenY);
    double sinvTem;
    sin>>sinvTem;
    this->setInvTem(sinvTem);
    int pNumUdv, pUdvIntv;
    sin>>pNumUdv>>pUdvIntv;
    this->setUdv(pNumUdv,pUdvIntv);
    int pSeed;
    sin>>pSeed;
    this->setSeed(pSeed);
    double pWin;
    sin>>pWin;
    this->setWin(pWin);
    double pHmltParm[NUM_HMLT_PARM];
    for(int i=0;i<NUM_HMLT_PARM;i++){
        sin>>pHmltParm[i];
    }
    this->model.setHmltParm(pHmltParm);
    for(int i=0;i<this->numSlc;i++){
        for(int j=0;j<this->model.numSite;j++){
            sin>>this->conf[i][j];
        }
    }
    sin.close();
    
    this->initMatK();
    this->prepScanFromZero();
}

void BscMC::calcNeqGrnKSngOrbt(MatrixXd* neqOrigR,double** neqK)const{
    int** spacRSiteList=new int*[this->model.numSite];
    for(int i=0;i<this->model.numSite;i++){
        spacRSiteList[i]=new int[this->model.numSite];
        this->findSpacRSite(i,spacRSiteList[i]);
    }

    double** neqR=new double*[this->numSlc];
    for(int i=0;i<this->numSlc;i++){
        neqR[i]=new double[this->model.numSite];
    }
    for(int slcInd=0;slcInd<this->numSlc;slcInd++){
        for(int rInd=0;rInd<this->model.numSite;rInd++){
            double sum=0;
            int j;
            for(int i=0;i<this->model.numSite;i++){
                j=spacRSiteList[rInd][i];
                sum+=neqOrigR[slcInd](i,j);
            }
            neqR[slcInd][rInd]=sum/this->model.numSite;
        }
    }

    double** cosList=new double*[this->model.numSite];
    for(int i=0;i<this->model.numSite;i++){
        cosList[i]=new double[this->model.numSite];
    }
    int kInd;
    int rInd;
    double dTmp;
    for(int kX=0;kX<this->model.lenX;kX++){
        for(int kY=0;kY<this->model.lenY;kY++){
            kInd=kX+kY*this->model.lenX;
            for(int rX=0;rX<this->model.lenX;rX++){
                for(int rY=0;rY<this->model.lenY;rY++){
                    rInd=rX+rY*this->model.lenX;
                    dTmp=double(kX*rX)/this->model.lenX+double(kY*rY)/this->model.lenY;
                    cosList[kInd][rInd]=cos(2*PI*dTmp);
                }
            }
        }
    }

    for(int kInd=0;kInd<this->model.numSite;kInd++){
        for(int slcInd=0;slcInd<this->numSlc;slcInd++){
            double sum=0;
            for(int rInd=0;rInd<this->model.numSite;rInd++){
                sum+=neqR[slcInd][rInd]*cosList[kInd][rInd];
            }
            neqK[kInd][slcInd]=sum;
        }
    }

    for(int i=0;i<this->model.numSite;i++){
        delete[] spacRSiteList[i];
    }
    delete[] spacRSiteList;
    for(int i=0;i<this->numSlc;i++){
        delete[] neqR[i];
    }
    delete[] neqR;
    for(int i=0;i<this->model.numSite;i++){
        delete[] cosList[i];
    }
    delete[] cosList;
}

void BscMC::calcNeqGrnK(MatrixXd* neqOrigRA,MatrixXd* neqOrigRB,double** neqKA,double** neqKB)const{
    this->calcNeqGrnKSngOrbt(neqOrigRA,neqKA);
    this->calcNeqGrnKSngOrbt(neqOrigRB,neqKB);
}

double BscMC::calcAvrgOcupSngSlc(MatrixXd& equGrnIn)const{
    return 1-equGrnIn.trace()/this->model.numSite;
}
double BscMC::calcAvrgOcup(MatrixXd* equGrnIn)const{
    double sum=0;
    for(int i=0;i<this->numSlc;i++){
        sum+=this->calcAvrgOcupSngSlc(equGrnIn[i]);
    }
    return sum/this->numSlc;
}

double BscMC::calcPairCorrS(MatrixXd* equGrnAIn,MatrixXd* equGrnBIn)const{
    int mX=this->model.lenX/2;
    int mY=this->model.lenY/2;
    double* sitePair=new double[this->model.numSite];
    for(int rX=0;rX<this->model.lenX;rX++){
        int iX=(rX+mX)%this->model.lenX;
        for(int rY=0;rY<this->model.lenY;rY++){
            int iY=(rY+mY)%this->model.lenY;
            sitePair[rX+rY*this->model.lenX]=iX+iY*this->model.lenX;            
        }
    }

    double sum=0;
    for(int slc=0;slc<this->numSlc;slc++){
        for(int site=0;site<this->model.numSite;site++){
            sum+=equGrnAIn[slc](site,sitePair[site]);
            sum+=equGrnBIn[slc](site,sitePair[site]);
        }
    }

    delete[] sitePair;
    double res=sum/this->numSlc/this->model.numSite/2;
    return res;
}

double BscMC::calcPDOS(double** neqKIn)const{
    double sum=0;
    int midSlc=this->numSlc/2;
    for(int kInd=0;kInd<this->model.numSite;kInd++){
        sum+=neqKIn[kInd][midSlc];
    }
    double res=sum*this->invTem/this->model.numSite;
    return res;
}

void BscMC::showOcupConf(const string filePath){
    ofstream dout;
    dout.open(filePath);

    MatrixXd* equTmpA=new MatrixXd[this->numSlc];
    MatrixXd* equTmpB=new MatrixXd[this->numSlc];

    this->mesrEquGrn(equTmpA,equTmpB);
    for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            dout<<1-equTmpA[0](siteInd,siteInd)<<" ";
        }
        for(int siteInd=0;siteInd<this->model.numSite;siteInd++){
            dout<<1-equTmpB[0](siteInd,siteInd)<<" ";
        }
        dout<<endl;


    delete[] equTmpA;
    delete[] equTmpB;
    dout.close();
}


#endif