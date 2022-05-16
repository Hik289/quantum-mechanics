#include "head.h"

void sphToCar(const double* const coorSph,double* coorCar){
    coorCar[0]=sin(coorSph[0])*cos(coorSph[1]);
    coorCar[1]=sin(coorSph[0])*sin(coorSph[1]);
    coorCar[2]=cos(coorSph[0]);
}

double calcACFun(const vector<double>& data,const int t){
    const int dataRan=data.size()-t;
    double res1=0;
    double res2=0;
    double res3=0;
    for(int i=0;i<dataRan;i++){
        res1+=data[i]*data[i+t];
        res2+=data[i];
        res3+=data[i+t];
    }
    res1=res1/dataRan;
    res2=res2/dataRan;
    res3=res3/dataRan;
    return res1-res2*res3;
}

double calcACTime(const vector<double>& data){
    vector<double> dataGood(0);
    double acFun;
    for(int t=0;t<data.size();t++){
        acFun=calcACFun(data,t);
        if(acFun>pow(0.1,9)){
            dataGood.push_back(acFun);
        }else{
            break;
        }
    }

    double resSum=0;
    for(int i=1;i<dataGood.size();i++){
        resSum+=(1-double(i)/dataGood.size())*dataGood[i];
    }
    double res;
    if(dataGood.size()==0){
        res=-1;
    }else{
        res=resSum/dataGood[0];
    }
    return res;
}


double calcMean(const vector<double>& data){
    double resSum=0;
    for(int i=0;i<data.size();i++){
        resSum+=data[i];
    }
    double res=resSum/data.size();
    return res;
}

double calcSqMean(const vector<double>& data){
    double sum=0;
    for(int i=0;i<data.size();i++){
        sum+=data[i]*data[i];
    }
    return sum/data.size();
}

double calcQuMean(const vector<double>& data){
    double sum=0;
    for(int i=0;i<data.size();i++){
        sum+=pow(data[i],4);
    }
    return sum/data.size();
}

double calcNaiErr(const vector<double>& data){
    const int len=data.size();
    double sumSq=0;
    for(int i=0;i<len;i++){
        sumSq+=data[i]*data[i];
    }
    const double meanData=calcMean(data);
    double resSq=(sumSq-len*meanData*meanData)/len/(len-1);
    return sqrt(resSq);
}

double calcErrBin(const vector<double>& data,const int numBin){
    const int lenBin=data.size()/numBin;
    vector<double> dataBin(numBin);
    for(int binInd=0;binInd<numBin;binInd++){
        double sumBin=0;
        for(int i=binInd*lenBin;i<(binInd+1)*lenBin;i++){
            sumBin+=data[i];
        }
        dataBin[binInd]=sumBin/lenBin;
    }
    return calcNaiErr(dataBin);
}

double calcErrAC(const vector<double>& data){
    const double acTime=calcACTime(data);
    const double naiErr=calcNaiErr(data);
    return naiErr*sqrt(2*acTime+1);
}

double calcVecNorm(const vector<double>& vec){
    double resSq=0;
    for(int i=0;i<vec.size();i++){
        resSq+=vec[i]*vec[i];
    }
    return sqrt(resSq);
}