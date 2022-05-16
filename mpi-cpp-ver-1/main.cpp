#include "head.h"

int main(){
    Tsk tsk;
    tsk.set(prmMinList,prmMaxList,numSlcList);
    int tagPrm=0;
    int tagStrFct=1;
    int tagEne=2;

    //MPI procedure starts here
    MPI_Init(NULL,NULL);
    int comm_size;
    int comm_rank;
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&comm_rank);

    if(comm_rank!=0){
        SysMC sysMC;
        int prmSetCode=-1;
        vector<double> prmSet(tsk.numPrm);
        double strFct[lat.numCell];
        int smltAnnlFlag=0;
        int resCode=-1;
        double resEne=0;

        while(1){
            MPI_Recv(&prmSetCode,1,MPI_INT,0,tagPrm,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            if(prmSetCode==-1){
                break;
            }
            tsk.codeToPrmSet(prmSetCode,prmSet);
            sysMC.setPrmHmlt(prmSet);
            sysMC.ranMT.setSeed(mtSeedGlo);
            sysMC.randomizeConf();
            sysMC.sprAng=PI;
            sysMC.smltAnnl(smltAnnlFlag,temIni,temEnd,sclFct,winAdjIntv);
            sysMC.calcStrFct(strFct);
            int judPeak=(strFct[6]>peakMin)&&(strFct[72]>peakMin)&&(strFct[78]>peakMin);
            // int judPeak=(strFct[6]>peakMin)||(strFct[72]>peakMin)||(strFct[78]>peakMin);
            int jud=judPeak*smltAnnlFlag;
            if(jud){
                resCode=prmSetCode;
                MPI_Send(&resCode,1,MPI_INT,0,tagPrm,MPI_COMM_WORLD);
                resEne=sysMC.eneTot();
                MPI_Send(&resEne,1,MPI_DOUBLE,0,tagEne,MPI_COMM_WORLD);
                MPI_Send(strFct,lat.numCell,MPI_DOUBLE,0,tagStrFct,MPI_COMM_WORLD);
            }else{
                resCode=-1;
                MPI_Send(&resCode,1,MPI_INT,0,tagPrm,MPI_COMM_WORLD);
            }
        }
    }else{
        clock_t startClock=clock();

        /*show procedure information*/
        ofstream iout;
        iout.open(outLoc+"info.dat");
        iout<<"================================================="<<endl;
        iout<<"lattice size:              "<<lat.numCellX<<" "<<lat.numCellY<<endl;
        iout<<"global random seed:        "<<mtSeedGlo<<endl;
        iout<<"number of slices:          ";
        for(int i=0;i<numSlcList.size();i++){
            iout<<numSlcList[i]<<" ";
        }
        iout<<endl;

        iout<<"left bound of parameters:  ";
        for(int i=0;i<prmMinList.size();i++){
            iout<<prmMinList[i]<<" ";
        }
        iout<<endl;

        iout<<"right bound of parameters: ";
        for(int i=0;i<prmMaxList.size();i++){
            iout<<prmMaxList[i]<<" ";
        }
        iout<<endl;
        iout<<"================================================="<<endl;
        iout<<endl;

        int prmSetCode=0;
        int prmSetCnt=0;
        int numPrmSet=tsk.calcNumPrmSet();
        int resCode;
        double resEne;
        double resStrFct[lat.numCell];
        MPI_Status stts;
        ofstream rout;
        rout.open(outLoc+"strFct.dat");

        /*publish initial tasks to all*/
        for(int thrInd=1;thrInd<comm_size;thrInd++){
            MPI_Send(&prmSetCode,1,MPI_INT,thrInd,tagPrm,MPI_COMM_WORLD);
            prmSetCnt++;
            if(prmSetCnt>=numPrmSet){
                prmSetCode=-1;
            }else{
                prmSetCode++;
            }
        }

        /*receive results and re-publish tasks*/
        while(prmSetCnt<(numPrmSet+comm_size-1)){
            MPI_Recv(&resCode,1,MPI_INT,MPI_ANY_SOURCE,tagPrm,MPI_COMM_WORLD,&stts);
            if(resCode>-1){
                MPI_Recv(&resEne,1,MPI_DOUBLE,stts.MPI_SOURCE,tagEne,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Recv(resStrFct,lat.numCell,MPI_DOUBLE,stts.MPI_SOURCE,tagStrFct,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                rout<<resCode<<" "<<resEne<<endl;
                vector<double> temp(tsk.numPrm);
                tsk.codeToPrmSet(resCode,temp);
                for(int i=0;i<7;i++){
                    rout<<temp[i]<<" ";
                }
                rout<<endl;
                for(int i=0;i<lat.numCell;i++){
                    rout<<resStrFct[i]<<" ";
                }
                rout<<endl;
                rout<<endl;
            }
            if( (prmSetCnt+1)%10000==0 ){
                iout<<(prmSetCnt+1)<<" parameters set calculated"<<endl;
            }
            MPI_Send(&prmSetCode,1,MPI_INT,stts.MPI_SOURCE,tagPrm,MPI_COMM_WORLD);
            prmSetCnt++;
            if(prmSetCnt>=numPrmSet){
                prmSetCode=-1;
            }else{
                prmSetCode++;
            }
        }

        rout.close();       
        clock_t endClock=clock();
        iout<<"time spent: "<<double(endClock-startClock)/CLOCKS_PER_SEC<<" seconds"<<endl;
        iout.close();
    }

    MPI_Finalize();
    return 0;
}