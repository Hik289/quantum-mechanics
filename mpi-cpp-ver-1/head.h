/**********************************************************
 * Using simulated annealing Monte Carlo to find parameters, with which the ground state is of triple-Q order.
 * Interactions involved (temporary): Kitaev, Gamma, 1st, 2nd, 3rd nearest neighbor Heisenberg.
 * Author: Kai-Wei Sun
 * Update date: 2020-11-30
 * **********************************************************/

#ifndef _HEAD_H
#define _HEAD_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <string>
#include <ctime>
#include <iomanip>
#include <algorithm>
#include <mpi.h>

#define PI 3.141592653589793

using namespace std;

class Lattice;
class ranMT;
class SysMC;
class Tsk;

extern const int numCellXGlo;
extern const int numCellYGlo;
extern const Lattice lat;

extern double temIni;
extern double temEnd;
extern double sclFct;
extern int winAdjIntv;
extern int mntrIntv;
extern int mtSeedGlo;

extern const vector<int> numSlcList;
extern const vector<double> prmMinList;
extern const vector<double> prmMaxList;
extern const double peakMin;
extern const string inLoc;
extern const string outLoc;

/*function declaration*/
/*transform from spherical coordinates to cartesian coordinates*/
void sphToCar(const double* const coorSph,double* coorCar);

/*data analysis*/
/*calculate autocorrelation time calculation*/
double calcACFun(const vector<double>& data,const int t);
/*calculate auto-correlation time*/
double calcACTime(const vector<double>& data);
/*calculate mean value given data in vector<double> form*/
double calcMean(const vector<double>& data);
/*calculate squared average*/
double calcSqMean(const vector<double>& data);
/*calculate quartic average*/
double calcQuMean(const vector<double>& data);
/*calculate naive error of data*/
/*accurate error is naive error multiplied by square root of (2*acTime+1)*/
double calcNaiErr(const vector<double>& data);
/*calculate error using binning analysis*/
/*this method is much faster than calculating auto-correlation time*/
double calcErrBin(const vector<double>& data,const int numBin);
/*calculate accurate error using auto-correlation time*/
/*thousand times slower than binning analysis*/
double calcErrAC(const vector<double>& data);
/*calculate 2-norm of vector<double>*/
double calcVecNorm(const vector<double>& vec);

class Lattice{
    public:
    int numCellX;
    int numCellY;
    int numCell;
    int numSite;

    int** bondList1X;
    int** bondList1Y;
    int** bondList1Z;
    int** bondList1;

    int** bondList2X;
    int** bondList2Y;
    int** bondList2Z;
    int** bondList2;

    int** bondList3X;
    int** bondList3Y;
    int** bondList3Z;
    int** bondList3;

    int** siteConnList1;
    int** siteConnList2;
    int** siteConnList3;

    /*constructor, copy constructor, overload "=" and destructor*/
    Lattice(const int numCellXIn=6,const int numCellYIn=6);
    Lattice(const Lattice& obj);
    Lattice& operator=(const Lattice& obj);
    ~Lattice();

    /*encode and decode*/
    void coorToCode(const int* const coor,int& code);
    void codeToCoor(const int code,int* coor);
};

class RanMT{
    public:
    RanMT(const unsigned int mtSeed=0);
    ~RanMT();

    /*read and reset seed*/
    unsigned int readSeed()const;
    void setSeed(const unsigned int mtSeed=0);
    /*return random integer*/
    int ranInt(const int ranMin=0,const int ranMax=1);
    /*return random real between ranMin to ranMax with default setting 0~1*/
    double ranReal(const double ranMin=0,const double ranMax=1);
    

    private:
    mt19937 mt;
    unsigned int mtMax;
    unsigned int curSeed;
};

class SysMC{
    public:
    double** conf;
    double invTem;
    double sprAng;
    RanMT ranMT;

    /*Hamiltonian parameters*/
    /*no minus sign in front of parameters in the Hamiltonian for convenience*/
    /*for example, J1<0 means FM*/
    double prmHsb1;
    double prmHsb2;
    double prmHsb3;
    double prmKtv;
    double prmGmm;
    /*asymmetric Gamma: s_i^\alpha s_j^\beta-s_i^\beta s_j^\alpha*/
    /*i belongs to sublattice A, j belongs to sublatticeB*/
    /*\gamma \alpha \beta can be xyz, yzx or zxy*/
    double prmGmmAsym;
    double prmAnIso;

    SysMC();
    SysMC(const SysMC& obj);
    SysMC& operator=(const SysMC& obj);
    ~SysMC();

    /*show parameters of Monte Carlo*/
    void showPrmMC(ostream& dest=cout)const;
    /*show parameters of Hamiltonian*/
    void showPrmHmlt(ostream& dest=cout)const;
    /*set parameters of Monte Carlo*/
    void setPrmMC(const double invTemIn,const double sprAngIn,const unsigned mtSeedIn);
    /*set parameters of Hamiltonian*/
    void setPrmHmlt(const vector<double>& prmHmlt,ostream& dest=cout);

    /*show configuration*/
    void showConf(ostream& dest=cout)const;
    void showConfCar(ostream& dest=cout)const;
    /*generate a random site*/
    int ranSite();
    /*generate random spin uniformly on the whole sphere*/
    void ranStateSph(double* state);
    /*generate a random spin uniformly distributed inside a cone given center spin and spread angle*/
    void ranStateCone(double* state,const double polCen,const double aziCen,const double sprAng);
    /*set configuration*/
    void setConf(double** confIn);
    /*randomize configuration*/
    void randomizeConf();

    /*total energy calculation*/
    double eneTotHsb1(double** confCar);
    double eneTotHsb2(double** confCar);
    double eneTotHsb3(double** confCar);
    double eneTotKtv(double** confCar);
    double eneTotGmm(double** confCar);
    double eneTotGmmAsym(double** confCar);
    double eneTotAnIso(double** confCar);
    double eneTot();
    /*calculate energy difference*/
    double eneDiff(const int siteUpd,double* spinUpd);

    /*try to update a spin*/
    void tryUpdSpin(const int siteUpd,double* spinUpd);
    /*overload: return accept flag*/
    void tryUpdSpin(const int siteUpd,double* spinUpd,int& accFlag);
    void tryUpdSpin(int& accFlag);
    /*sweep in two manners: sequential- and random-site-chosen*/
    void swpRan();
    void swpSqn();
    /*overload: calculate accept rate*/
    void swpRan(double& accRate);
    void swpSqn(double& accRate);
    /*for these two ways, calculate accept ratio*/
    double calcAccRateFast();
    double calcAccRateRan();
    double calcAccRateSqn();
    /*adjust spread angle*/
    void adjustWinFast();
    void adjustWinRan(ostream& iout=cout);
    void adjustWinSqn(ostream& iout=cout);

    /*calclate static structure factor*/
    /*each line: kx,ky,strFct*/
    void calcStrFct(double** res);
    /*only strFct*/
    void calcStrFct(double* res);
    /*overload: calculate structure factor and export*/
    void calcStrFct(ostream& iout=cout);

    /*simulated annealing*/
    /*overload: naive version of simulated annealing*/
    void smltAnnl(vector<double>& invTemPath,const int numSwp,ostream& iout=cout);
    /*overload: thermodynamic simulated annealing*/
    void smltAnnl(const double temIni=10000,const double temEnd=0.001,const double sclFct=10,const int winAdjIntv=10000,const int mntrIntv=10000,ostream& iout=cout);
    /*overload: thermodynamic simulated annealing with no output info*/
    /*if not found, return flag: succeed=1, failed=0*/
    void smltAnnl(int& flag,const double temIni=10000,const double temEnd=0.001,const double sclFct=10,const int winAdjIntv=10000);
};

class Tsk{
    public:
    int numPrm;
    double** prmStng;
    int* numSlcList;

    Tsk(const int numPrmIn=7);
    Tsk(const Tsk& obj);
    Tsk& operator=(const Tsk& obj);
    ~Tsk();

    void set(const vector<double>& prmMinList,const vector<double>& prmMaxList,const vector<int>& numVarList);
    void codeToPrmSet(const int code,vector<double>& prm)const;
    int calcNumPrmSet()const;
};



#endif