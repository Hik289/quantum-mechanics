/*
Basic data analyzing functions.
*/

#ifndef _DATA_ANALYZER_HPP
#define _DATA_ANALYZER_HPP

#include <vector>
#include <cmath>
using std::vector;

/*calculate average of given data*/
double calcAvrg(const vector<double>& data);
/*calculate squared average*/
double calcSqAvrg(const vector<double>& data);
/*calculate quartic average*/
double calcQtAvrg(const vector<double>& data);
/*calculate auto-correlation function*/
double calcAutoCorrFunc(const vector<double>& data,const int t);
/*calculate auto-correlation time*/
double calcAutoCorrTime(const vector<double>& data);
/*calculate naive error ignoring auto-correlation*/
double calcNaiveErr(const vector<double>& data);
/*calculate error in the binning manner*/
double calcErrBin(const vector<double>& data,const int numBin=5);
/*calculate error by calculating auto-correlation time*/
double calcErrAC(const vector<double>& data);
/*clculate norm of a given vector*/
double calcVecNorm(const vector<double>& vec);

double calcAvrg(const vector<double>& data){
    double resSum=0;
    for(int i=0;i<data.size();i++){
        resSum+=data[i];
    }
    double res=resSum/data.size();
    return res;
}

double calcSqAvrg(const vector<double>& data){
    double sum=0;
    for(int i=0;i<data.size();i++){
        sum+=data[i]*data[i];
    }
    return sum/data.size();
}

double calcQtAvrg(const vector<double>& data){
    double sum=0;
    for(int i=0;i<data.size();i++){
        sum+=pow(data[i],4);
    }
    return sum/data.size();
}

double calcAutoCorrFunc(const vector<double>& data,const int t){
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

double calcAutoCorrTime(const vector<double>& data){
    vector<double> dataGood(0);
    double acFun;
    for(int t=0;t<data.size();t++){
        acFun=calcAutoCorrFunc(data,t);
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

double calcNaiveErr(const vector<double>& data){
    const int len=data.size();
    double sumSq=0;
    for(int i=0;i<len;i++){
        sumSq+=data[i]*data[i];
    }
    const double meanData=calcAvrg(data);
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
    return calcNaiveErr(dataBin);
}

double calcErrAC(const vector<double>& data){
    const double acTime=calcAutoCorrTime(data);
    const double naiErr=calcNaiveErr(data);
    return naiErr*sqrt(2*acTime+1);
}

double calcVecNorm(const vector<double>& vec){
    double resSq=0;
    for(int i=0;i<vec.size();i++){
        resSq+=vec[i]*vec[i];
    }
    return sqrt(resSq);
}

#endif