/*initialization: parameter setting and equilibrium*/
/*the path is recorded to ensure equilibrium*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <mpi.h>
#include "model.hpp"
#include "data_analyzer.hpp"
#include "random_Mersenne_twister.hpp"
#include "Eigen/Eigen"
#include "bscMC.hpp"

using namespace std;
using namespace Eigen;

/*gobal parameters*/
string dirData;
int tagFinish=0;
int tagTaskInvTem=1;
int tagTaskCpl=2;
int tagTaskSeed=3;

/*auxiliary functions declaration*/
void parseCmd(int argc, char** argv);
void readParm(BscMC& obj,const string fileParm);

int main(int argc,char** argv){
    MPI_Init(&argc,&argv);
    int comm_size;
    int comm_rank;
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&comm_rank);

    parseCmd(argc,argv);
    stringstream ss;
    string strRank;
    ss<<comm_rank;
    ss>>strRank;
    ss.clear();

    if(comm_rank==0){
        ofstream lout;
        lout.open(dirData+"/log/log_"+strRank+".dat");
        clock_t clk1=clock();

        /*read tasks and publish to workers*/
        ifstream tin;
        tin.open(dirData+"/task.dat");
        double invTemTmp;
        double cplTmp;
        int seedTmp;
        for(int thd=1;thd<comm_size;thd++){
            tin>>invTemTmp;
            MPI_Send(&invTemTmp,1,MPI_DOUBLE,thd,tagTaskInvTem,MPI_COMM_WORLD);
            tin>>cplTmp;
            MPI_Send(&cplTmp,1,MPI_DOUBLE,thd,tagTaskCpl,MPI_COMM_WORLD);
            tin>>seedTmp;
            MPI_Send(&seedTmp,1,MPI_INT,thd,tagTaskSeed,MPI_COMM_WORLD);
        }
        tin.close();

        /*receive task-finish-signals*/
        MPI_Status mpiStts;
        int flagFinish;
        for(int thd=1;thd<comm_size;thd++){
            MPI_Recv(&flagFinish,1,MPI_INT,MPI_ANY_SOURCE,tagFinish,MPI_COMM_WORLD,&mpiStts);
            lout<<mpiStts.MPI_SOURCE<<" finished"<<endl;
        }
        clock_t clk2=clock();
        lout<<"time spent: "<<double(clk2-clk1)/CLOCKS_PER_SEC<<" seconds"<<endl;
        lout.close();
    }else{
        ofstream lout;
        lout.open(dirData+"/log/log_"+strRank+".dat");

        BscMC bscMC;
        readParm(bscMC,dirData+"/parm.dat");
        double invTemTmp;
        double cplTmp;
        int seedTmp;
        MPI_Recv(&invTemTmp,1,MPI_DOUBLE,0,tagTaskInvTem,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(&cplTmp,1,MPI_DOUBLE,0,tagTaskCpl,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(&seedTmp,1,MPI_INT,0,tagTaskSeed,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        bscMC.setInvTem(invTemTmp);
        bscMC.model.strCpl=cplTmp;
        bscMC.setSeed(seedTmp);
        bscMC.ranConf();
        bscMC.showParm(lout);
        lout.close();

        bscMC.exportStts(dirData+"/stts/stts_"+strRank+".dat");

        int flagFinish=1;
        MPI_Send(&flagFinish,1,MPI_INT,0,tagFinish,MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}



void parseCmd(int argc, char** argv){
    int status=0;
    for(int i=1;i<argc;i++){
        switch(status){
            case 0:
            if(argv[i][0]=='-'){
                switch(argv[i][1]){
                    case 'd': 
                    status=1;
                    break;
                    default:
                    cout<<"ERROR: command not recognized."<<endl;
                    return;
                }
            }
            case 1:
            dirData=argv[i];
            status=0;
        }
    }
}

void readParm(BscMC& obj,const string fileParm){
    ifstream pin;
    pin.open(fileParm);
    int pLenX,pLenY;
    pin>>pLenX>>pLenY;
    obj.setLattSize(pLenX,pLenY);
    double pInvTem;
    pin>>pInvTem;
    obj.setInvTem(pInvTem);
    int pNumUdv, pUdvIntv;
    pin>>pNumUdv>>pUdvIntv;
    obj.setUdv(pNumUdv,pUdvIntv);
    int pSeed;
    pin>>pSeed;
    obj.setSeed(pSeed);
    double pWin;
    pin>>pWin;
    obj.setWin(pWin);
    double pHmltParm[NUM_HMLT_PARM];
    for(int i=0;i<NUM_HMLT_PARM;i++){
        pin>>pHmltParm[i];
    }
    obj.model.setHmltParm(pHmltParm);
    pin.close();
}