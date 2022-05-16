#include "head.h"

int main(){
    SysMC sysMC;

    // this->prmHsb1=prmHmlt[0];
    // this->prmHsb2=prmHmlt[1];
    // this->prmHsb3=prmHmlt[2];
    // this->prmKtv=prmHmlt[3];
    // this->prmGmm=prmHmlt[4];
    // this->prmGmmAsym=prmHmlt[5];
    // this->prmAnIso=prmHmlt[6];

    // vector<double> prmHmlt={0,0,0,1,-1,0,0};
    vector<double> prmHmlt={-1.28205,-0.769231,1.38462,-1,2,0,0};
    sysMC.setPrmHmlt(prmHmlt);
    sysMC.randomizeConf();



    double temIni=10;
    double temEnd=0.0001;
    double sclFct=1.01;
    int winAdjIntv=10000;
    int mntrIntv=100000;

    sysMC.ranMT.setSeed(0);
    sysMC.sprAng=1;
    sysMC.smltAnnl(temIni,temEnd,sclFct,winAdjIntv,mntrIntv);

    ofstream fout;
    fout.open("/home/kevinsun/work/projects/tripleQ/cpp-ver-1/dataOut/conf.dat");
    sysMC.showConfCar(fout);
    fout.close();
    fout.open("/home/kevinsun/work/projects/tripleQ/cpp-ver-1/dataOut/strFct.dat");
    sysMC.calcStrFct(fout);
    fout.close();

    return 0;
}