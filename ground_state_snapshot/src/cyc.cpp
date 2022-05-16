/*measurement*/
/*configuration and parameters to be imported from stts_*.dat*/

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
int numCyc;
int tagFinish=0;

/*auxiliary functions declaration*/
void parseCmd(int argc, char** argv);
void readParm(const string fileParm);

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
        lout.open(dirData+"/log/log_"+strRank+".dat",ios::app);
        clock_t clk1=clock();
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
        lout.open(dirData+"/log/log_"+strRank+".dat",ios::app);
        lout<<"----------------------------------------------"<<endl;

        BscMC bscMC;
        bscMC.importStts(dirData+"/stts/stts_"+strRank+".dat");
        readParm(dirData+"/parm.dat");

        bscMC.printRoundOffInfo(lout);

        int rprtIntv=numCyc/10;
        int winAdjIntv=50;
        int mesrIntv=50;
        MatrixXd* neqRawL=new MatrixXd[bscMC.numSlc];
        MatrixXd* neqRawS=new MatrixXd[bscMC.numSlc];

        MatrixXd* neqRawAL=new MatrixXd[bscMC.numSlc];
        MatrixXd* neqRawAS=new MatrixXd[bscMC.numSlc];
        MatrixXd* neqRawBL=new MatrixXd[bscMC.numSlc];
        MatrixXd* neqRawBS=new MatrixXd[bscMC.numSlc];
        for(int i=0;i<bscMC.numSlc;i++){
            neqRawL[i]=MatrixXd(bscMC.model.numSite,bscMC.model.numSite);
            neqRawS[i]=MatrixXd(bscMC.model.numSite,bscMC.model.numSite);
            neqRawL[i].setZero();
            neqRawS[i].setZero();

            neqRawAL[i]=MatrixXd(bscMC.model.numSite,bscMC.model.numSite);
            neqRawAS[i]=MatrixXd(bscMC.model.numSite,bscMC.model.numSite);
            neqRawBL[i]=MatrixXd(bscMC.model.numSite,bscMC.model.numSite);
            neqRawBS[i]=MatrixXd(bscMC.model.numSite,bscMC.model.numSite);
            neqRawAL[i].setZero();
            neqRawAS[i].setZero();
            neqRawBL[i].setZero();
            neqRawBS[i].setZero();
        }
        MatrixXd* equTmpA=new MatrixXd[bscMC.numSlc];
        MatrixXd* equTmpB=new MatrixXd[bscMC.numSlc];
        MatrixXd* neqTmpA=new MatrixXd[bscMC.numSlc];
        MatrixXd* neqTmpB=new MatrixXd[bscMC.numSlc];

        ofstream eout;
        eout.open(dirData+"/equ_check/path_"+strRank+".dat",ios::app);
        int numMesr=0;

        int numAL=0;
        int numBL=0;
        double sumOcupL=0;
        double sumOcupS=0;
        double sumOrdFM=0;
        double sumCorrFM=0;
        double sumFld=0;
        double sumFldSq=0;
        double sumFldQt=0;
        double sumPairCorrS=0;
        for(int cycInd=0;cycInd<numCyc;cycInd++){
            bscMC.cycScan();
            if( (cycInd+1)%winAdjIntv==0 ){
                bscMC.adjWin();
            }
            if( (cycInd+1)%rprtIntv==0 ){
                lout<<cycInd+1<<" out of "<<numCyc<<" steps finished"<<endl;
            }
            if( (cycInd+1)%mesrIntv==0 ){
                numMesr++;
                /*equal-time*/
                bscMC.mesrEquGrn(equTmpA,equTmpB);
                double ocupA=bscMC.calcAvrgOcup(equTmpA);
                double ocupB=bscMC.calcAvrgOcup(equTmpB);
                double ordFM=bscMC.calcOrdFM(equTmpA,equTmpB);
                sumOrdFM+=abs(ordFM);
                sumCorrFM+=bscMC.calcEquCorrFM(equTmpA,equTmpB);
                sumFld+=abs(bscMC.mesrConjFld());
                sumFldSq+=bscMC.mesrConjFldSq();
                sumFldQt+=bscMC.mesrConjFldQt();
                sumPairCorrS+=bscMC.calcPairCorrS(equTmpA,equTmpB);

                /*non-equal time*/
                bscMC.mesrNeqGrn(neqTmpA,neqTmpB);

                if(ordFM>0){
                    numAL++;
                    sumOcupL+=ocupA;
                    sumOcupS+=ocupB;
                    for(int i=0;i<bscMC.numSlc;i++){
                        neqRawL[i]+=neqTmpA[i];
                        neqRawS[i]+=neqTmpB[i];

                        neqRawAL[i]+=neqTmpA[i];
                        neqRawBS[i]+=neqTmpB[i];
                    }
                }else{
                    numBL++;
                    sumOcupL+=ocupB;
                    sumOcupS+=ocupA;
                    for(int i=0;i<bscMC.numSlc;i++){
                        neqRawS[i]+=neqTmpA[i];
                        neqRawL[i]+=neqTmpB[i];

                        neqRawAS[i]+=neqTmpA[i];
                        neqRawBL[i]+=neqTmpB[i];
                    }
                }
                eout<<numMesr<<" "<<numAL<<" "<<numBL<<" ";
                eout<<ocupA<<" "<<ocupB<<" "<<ordFM<<endl;
            }
        }
        eout.close();

        for(int i=0;i<bscMC.numSlc;i++){
            neqRawL[i]=neqRawL[i]/numMesr;
            neqRawS[i]=neqRawS[i]/numMesr;

            neqRawAL[i]=neqRawAL[i]/numAL;
            neqRawBS[i]=neqRawBS[i]/numAL;
            neqRawAS[i]=neqRawAS[i]/numBL;
            neqRawBL[i]=neqRawBL[i]/numBL;
        }
        double** neqKL=new double*[bscMC.model.numSite];
        double** neqKS=new double*[bscMC.model.numSite];

        double** neqKAL=new double*[bscMC.model.numSite];
        double** neqKAS=new double*[bscMC.model.numSite];
        double** neqKBL=new double*[bscMC.model.numSite];
        double** neqKBS=new double*[bscMC.model.numSite];
        for(int i=0;i<bscMC.model.numSite;i++){
            neqKL[i]=new double[bscMC.numSlc];
            neqKS[i]=new double[bscMC.numSlc];

            neqKAL[i]=new double[bscMC.numSlc];
            neqKAS[i]=new double[bscMC.numSlc];
            neqKBL[i]=new double[bscMC.numSlc];
            neqKBS[i]=new double[bscMC.numSlc];
        }
        bscMC.calcNeqGrnK(neqRawL,neqRawS,neqKL,neqKS);

        bscMC.calcNeqGrnK(neqRawAL,neqRawAS,neqKAL,neqKAS);
        bscMC.calcNeqGrnK(neqRawBL,neqRawBS,neqKBL,neqKBS);

        ofstream rout;
        rout.open(dirData+"/res/res_"+strRank+".dat",ios::app);
        rout<<sumOcupL/numMesr<<" ";
        rout<<sumOcupS/numMesr<<" ";
        rout<<sumOrdFM/numMesr<<" ";
        rout<<sumCorrFM/numMesr<<" ";
        rout<<sumFld/numMesr<<" ";
        rout<<sumFldSq/numMesr<<" ";
        rout<<sumFldQt/numMesr<<" ";
        rout<<sumPairCorrS/numMesr<<endl;
        rout.close();
        rout.clear();

        rout.open(dirData+"/neq/neq_"+strRank+"_L.dat",ios::app);
        rout<<numMesr<<endl;
        for(int kInd=0;kInd<bscMC.model.numSite;kInd++){
            for(int slcInd=0;slcInd<bscMC.numSlc;slcInd++){
                rout<<neqKL[kInd][slcInd]<<" ";
            }
            rout<<endl;
        }
        rout.close();
        rout.clear();
        rout.open(dirData+"/neq/neq_"+strRank+"_S.dat",ios::app);
        rout<<numMesr<<endl;
        for(int kInd=0;kInd<bscMC.model.numSite;kInd++){
            for(int slcInd=0;slcInd<bscMC.numSlc;slcInd++){
                rout<<neqKS[kInd][slcInd]<<" ";
            }
            rout<<endl;
        }
        rout.close();
        rout.clear();

        rout.open(dirData+"/neq/neq_"+strRank+"_AL.dat",ios::app);
        rout<<numAL<<endl;
        for(int kInd=0;kInd<bscMC.model.numSite;kInd++){
            for(int slcInd=0;slcInd<bscMC.numSlc;slcInd++){
                rout<<neqKAL[kInd][slcInd]<<" ";
            }
            rout<<endl;
        }
        rout.close();
        rout.clear();
        rout.open(dirData+"/neq/neq_"+strRank+"_BS.dat",ios::app);
        rout<<numAL<<endl;
        for(int kInd=0;kInd<bscMC.model.numSite;kInd++){
            for(int slcInd=0;slcInd<bscMC.numSlc;slcInd++){
                rout<<neqKBS[kInd][slcInd]<<" ";
            }
            rout<<endl;
        }
        rout.close();
        rout.clear();
        rout.open(dirData+"/neq/neq_"+strRank+"_BL.dat",ios::app);
        rout<<numBL<<endl;
        for(int kInd=0;kInd<bscMC.model.numSite;kInd++){
            for(int slcInd=0;slcInd<bscMC.numSlc;slcInd++){
                rout<<neqKBL[kInd][slcInd]<<" ";
            }
            rout<<endl;
        }
        rout.close();
        rout.clear();
        rout.open(dirData+"/neq/neq_"+strRank+"_AS.dat",ios::app);
        rout<<numBL<<endl;
        for(int kInd=0;kInd<bscMC.model.numSite;kInd++){
            for(int slcInd=0;slcInd<bscMC.numSlc;slcInd++){
                rout<<neqKAS[kInd][slcInd]<<" ";
            }
            rout<<endl;
        }
        rout.close();

        delete[] neqRawL;
        delete[] neqRawS;

        delete[] neqRawAL;
        delete[] neqRawAS;
        delete[] neqRawBL;
        delete[] neqRawBS;
        delete[] neqTmpA;
        delete[] neqTmpB;
        delete[] equTmpA;
        delete[] equTmpB;
        for(int i=0;i<bscMC.model.numSite;i++){
            delete[] neqKL[i];
            delete[] neqKS[i];

            delete[] neqKAL[i];
            delete[] neqKAS[i];
            delete[] neqKBL[i];
            delete[] neqKBS[i];
        }
        delete[] neqKL;
        delete[] neqKS;

        delete[] neqKAL;
        delete[] neqKAS;
        delete[] neqKBL;
        delete[] neqKBS;
        lout.close();


        bscMC.showOcupConf(dirData+"/log/ocup_"+strRank+".dat");
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

void readParm(const string fileParm){
    ifstream pin;
    pin.open(fileParm);
    double dTmp;
    for(int i=0;i<16;i++){
        pin>>dTmp;
    }
    pin>>numCyc;
    pin.close();
}