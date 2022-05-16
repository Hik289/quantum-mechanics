#include "head.h"

/*lattice parameters*/
const int numCellXGlo=12;
const int numCellYGlo=12;
const Lattice lat(numCellXGlo,numCellYGlo);

/*parameters for simulated annealing*/
double temIni=10;
double temEnd=0.001;
double sclFct=10;
int winAdjIntv=10000;
int mtSeedGlo=0;

/*parameters for parameter searching*/
const vector<int> numSlcList={40,40,40,1,40,1,1};
const vector<double> prmMinList={-2,-2,-2,-1,-2,0,0};
const vector<double> prmMaxList={2,2,2,-1,2,0,0};
const double peakMin=-1;

const string inLoc="/home/kevinsun/work/projects/tripleQ/mpi-cpp-ver-1/dataIn/";
const string outLoc="/home/kevinsun/work/projects/tripleQ/mpi-cpp-ver-1/dataOut/";