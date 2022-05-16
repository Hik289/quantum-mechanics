/*parameters setting for nematicity model*/
#ifndef _MODEL_HPP
#define _MODEL_HPP

#include <iostream>
#include <fstream>

#define NUM_HMLT_PARM 9

using std::ostream;
using std::cout;
using std::endl;

class Model{
public:
    /*Hamiltonian parameters*/
    double strHop;
    double ratHop;
    double strCpl;
    double strDyn;
    double strGrad;
    double strSq;
    double strQt;
    double chem;
    double biasChem;
    /*lattice parameters*/
    int lenX;
    int lenY;
    int numSite;
    /*useful lattice lists*/
    int** bondX;
    int** bondY;
    int** fourSite1;
    int** fourSite2;

    /*construct and destruct*/
    Model(const int lenXIn=4,const int lenYIn=4);
    Model(const Model& obj);
    Model& operator = (const Model& obj);
    ~Model();

    /*parameter setting and reading*/
    void setHmltParm(double* parm);
    void setLattSize(const int lenXIn,const int lenYIn);
    void showModelParm(ostream& iout=cout)const;
    void showBond(ostream& iout=cout)const;
    void showFoutSite(ostream& iout=cout)const;


    /*auxiliary functions*/
    void coorToCode(const int* coorIn,int &codeOut)const;
    void codeToCoor(const int codeIn,int* coorOut)const;
    void genNeigSite(const int codeCen,int* codeNeig)const;
};

Model::Model(const int lenXIn,const int lenYIn){
    this->strHop=0;
    this->ratHop=0;
    this->strCpl=0;
    this->strDyn=0;
    this->strGrad=0;
    this->strSq=0;
    this->strQt=0;
    this->chem=0;
    this->biasChem=0;
    
    this->lenX=lenXIn;
    this->lenY=lenYIn;
    this->numSite=lenXIn*lenYIn;

    /*construct bondX and bondY*/
    this->bondX=new int*[this->numSite];
    this->bondY=new int*[this->numSite];
    for(int i=0;i<this->numSite;i++){
        this->bondX[i]=new int[2];
        this->bondY[i]=new int[2];
    }

    int coorTemp[2];
    int codeTemp;
    int bondInd;

    for(int xInd=0;xInd<this->lenX;xInd++){
        for(int yInd=0;yInd<this->lenY;yInd++){
            bondInd=yInd*this->lenX+xInd;
            coorTemp[0]=xInd;
            coorTemp[1]=yInd;
            this->coorToCode(coorTemp,codeTemp);
            this->bondX[bondInd][0]=codeTemp;
            this->bondY[bondInd][0]=codeTemp;

            coorTemp[0]=(xInd+1)%this->lenX;
            coorTemp[1]=yInd;
            this->coorToCode(coorTemp,codeTemp);
            this->bondX[bondInd][1]=codeTemp;

            coorTemp[0]=xInd;
            coorTemp[1]=(yInd+1)%this->lenY;
            this->coorToCode(coorTemp,codeTemp);
            this->bondY[bondInd][1]=codeTemp;
        }
    }

    /*construct fourSite1 and fourSite2*/
    this->fourSite1=new int*[this->numSite/4];
    this->fourSite2=new int*[this->numSite/4];
    for(int i=0;i<this->numSite/4;i++){
        this->fourSite1[i]=new int[4];
        this->fourSite2[i]=new int[4];
    }
    for(int xInd=0;xInd<this->lenX;xInd+=2){
        for(int yInd=0;yInd<this->lenY;yInd+=2){
            int coor0[2]={xInd,yInd};
            int coor1[2]={(xInd+1)%this->lenX,yInd};
            int coor2[2]={(xInd+1)%this->lenX,(yInd+1)%this->lenY};
            int coor3[2]={xInd,(yInd+1)%this->lenY};
            int code0,code1,code2,code3;
            this->coorToCode(coor0,code0);
            this->coorToCode(coor1,code1);
            this->coorToCode(coor2,code2);
            this->coorToCode(coor3,code3);
            this->fourSite1[xInd/2+yInd*this->lenX/4][0]=code0;
            this->fourSite1[xInd/2+yInd*this->lenX/4][1]=code1;
            this->fourSite1[xInd/2+yInd*this->lenX/4][2]=code2;
            this->fourSite1[xInd/2+yInd*this->lenX/4][3]=code3;
        }
    }
    for(int xInd=1;xInd<this->lenX;xInd+=2){
        for(int yInd=1;yInd<this->lenY;yInd+=2){
            int coor0[2]={xInd,yInd};
            int coor1[2]={(xInd+1)%lenX,yInd};
            int coor2[2]={(xInd+1)%lenX,(yInd+1)%this->lenY};
            int coor3[2]={xInd,(yInd+1)%this->lenY};
            int code0,code1,code2,code3;
            this->coorToCode(coor0,code0);
            this->coorToCode(coor1,code1);
            this->coorToCode(coor2,code2);
            this->coorToCode(coor3,code3);
            this->fourSite2[(xInd-1)/2+(yInd-1)*this->lenX/4][0]=code0;
            this->fourSite2[(xInd-1)/2+(yInd-1)*this->lenX/4][1]=code1;
            this->fourSite2[(xInd-1)/2+(yInd-1)*this->lenX/4][2]=code2;
            this->fourSite2[(xInd-1)/2+(yInd-1)*this->lenX/4][3]=code3;
        }
    }
}

Model::Model(const Model& obj){
    this->strHop=obj.strHop;
    this->ratHop=obj.ratHop;
    this->strCpl=obj.strCpl;
    this->strDyn=obj.strDyn;
    this->strGrad=obj.strGrad;
    this->strSq=obj.strSq;
    this->strQt=obj.strQt;
    this->chem=obj.chem;
    this->biasChem=obj.biasChem;
    this->lenX=obj.lenX;
    this->lenY=obj.lenY;
    this->numSite=this->lenX*this->lenY;
    
    this->bondX=new int*[this->numSite];
    this->bondY=new int*[this->numSite];
    for(int i=0;i<numSite;i++){
        this->bondX[i]=new int[2];
        this->bondY[i]=new int[2];
    }
    this->fourSite1=new int*[this->numSite/4];
    this->fourSite2=new int*[this->numSite/4];
    for(int i=0;i<this->numSite/4;i++){
        this->fourSite1[i]=new int[4];
        this->fourSite2[i]=new int[4];
    }

    for(int i=0;i<this->numSite;i++){
        for(int j=0;j<2;j++){
            this->bondX[i][j]=obj.bondX[i][j];
            this->bondY[i][j]=obj.bondY[i][j];
        }
    }

    for(int i=0;i<(numSite/4);i++){
        for(int j=0;j<4;j++){
            this->fourSite1[i][j]=obj.fourSite1[i][j];
            this->fourSite2[i][j]=obj.fourSite2[i][j];
        }
    }
}

Model& Model::operator= (const Model& obj){
    if(this==&obj){
        return *this;
    }
    /*here we need to deallocate bondX, bondY, fourSite1 and fourSite2*/
    for(int i=0;i<this->numSite;i++){
        delete[] bondX[i];
        delete[] bondY[i];
    }
    delete[] bondX;
    delete[] bondY;
    for(int i=0;i<this->numSite/4;i++){
        delete[] fourSite1[i];
        delete[] fourSite2[i];
    }
    delete[] fourSite1;
    delete[] fourSite2;

    /*assign new values*/
    this->strHop=obj.strHop;
    this->ratHop=obj.ratHop;
    this->strCpl=obj.strCpl;
    this->strDyn=obj.strDyn;
    this->strGrad=obj.strGrad;
    this->strSq=obj.strSq;
    this->strQt=obj.strQt;
    this->chem=obj.chem;
    this->biasChem=obj.biasChem;
    this->lenX=obj.lenX;
    this->lenY=obj.lenY;
    this->numSite=this->lenX*this->lenY;
    
    this->bondX=new int*[this->numSite];
    this->bondY=new int*[this->numSite];
    for(int i=0;i<numSite;i++){
        this->bondX[i]=new int[2];
        this->bondY[i]=new int[2];
    }
    this->fourSite1=new int*[this->numSite/4];
    this->fourSite2=new int*[this->numSite/4];
    for(int i=0;i<this->numSite/4;i++){
        this->fourSite1[i]=new int[4];
        this->fourSite2[i]=new int[4];
    }

    for(int i=0;i<this->numSite;i++){
        for(int j=0;j<2;j++){
            this->bondX[i][j]=obj.bondX[i][j];
            this->bondY[i][j]=obj.bondY[i][j];
        }
    }

    for(int i=0;i<(numSite/4);i++){
        for(int j=0;j<4;j++){
            this->fourSite1[i][j]=obj.fourSite1[i][j];
            this->fourSite2[i][j]=obj.fourSite2[i][j];
        }
    }

    return *this;
}

Model::~Model(){
    /*here we need to deallocate bondX and bondY*/
    for(int i=0;i<this->numSite;i++){
        delete[] bondX[i];
        delete[] bondY[i];
    }
    delete[] bondX;
    delete[] bondY;
    /*deallocate fourSite1 and fourSite2*/
    for(int i=0;i<this->numSite/4;i++){
        delete[] fourSite1[i];
        delete[] fourSite2[i];
    }
    delete[] fourSite1;
    delete[] fourSite2;
}

/*mutually transformation between code and coordinates*/
void Model::coorToCode(const int* coorIn,int &codeOut)const{
    codeOut=coorIn[0]+coorIn[1]*this->lenX;
}
void Model::codeToCoor(const int codeIn,int* coorOut)const{
    int y=codeIn/this->lenX;
    int x=codeIn-y*this->lenX;
    coorOut[0]=x;
    coorOut[1]=y;
}

/*return codes of connected sites (nearest neighbors)*/
void Model::genNeigSite(const int codeCen,int* codeNeig)const{
    int coorCen[2];
    this->codeToCoor(codeCen,coorCen);

    int coorTemp[2];
    int codeTemp;

    coorTemp[0]=coorCen[0];
    coorTemp[1]=(coorCen[1]-1+this->lenY)%this->lenY;
    this->coorToCode(coorTemp,codeTemp);
    codeNeig[0]=codeTemp;
    coorTemp[1]=(coorCen[1]+1)%this->lenY;
    this->coorToCode(coorTemp,codeTemp);
    codeNeig[1]=codeTemp;

    coorTemp[0]=(coorCen[0]-1+this->lenX)%this->lenX;
    coorTemp[1]=coorCen[1];
    this->coorToCode(coorTemp,codeTemp);
    codeNeig[2]=codeTemp;
    coorTemp[0]=(coorCen[0]+1)%this->lenX;
    this->coorToCode(coorTemp,codeTemp);
    codeNeig[3]=codeTemp;
}

void Model::showBond(ostream& iout)const{
    iout<<"===================bondX==================="<<endl;
    for(int i=0;i<this->numSite;i++){
        iout<<i+1<<": "<<this->bondX[i][0]<<" "<<this->bondX[i][1]<<endl;
    }
    iout<<"===================bondY==================="<<endl;
    for(int i=0;i<this->numSite;i++){
        iout<<i+1<<": "<<this->bondY[i][0]<<" "<<this->bondY[i][1]<<endl;
    }
}

void Model::showFoutSite(ostream& iout)const{
    iout<<"===================fourSite1==================="<<endl;
    for(int i=0;i<this->numSite/4;i++){
        iout<<i+1<<": ";
        for(int j=0;j<4;j++){
            iout<<this->fourSite1[i][j]<<" ";
        }
        iout<<endl;
    }
    iout<<"===================fourSite2==================="<<endl;
    for(int i=0;i<this->numSite/4;i++){
        iout<<i+1<<": ";
        for(int j=0;j<4;j++){
            iout<<this->fourSite2[i][j]<<" ";
        }
        iout<<endl;
    }
}

void Model::setHmltParm(double* parm){
    this->strHop=parm[0];
    this->ratHop=parm[1];
    this->strCpl=parm[2];
    this->strDyn=parm[3];
    this->strGrad=parm[4];
    this->strSq=parm[5];
    this->strQt=parm[6];
    this->chem=parm[7];
    this->biasChem=parm[8];
}

void Model::setLattSize(const int lenXIn,const int lenYIn){
    if( (lenXIn==this->lenX)&&(lenYIn==this->lenY) ){
        return;
    }

    for(int i=0;i<this->numSite;i++){
        delete[] bondX[i];
        delete[] bondY[i];
    }
    delete[] bondX;
    delete[] bondY;
    for(int i=0;i<this->numSite/4;i++){
        delete[] fourSite1[i];
        delete[] fourSite2[i];
    }
    delete[] fourSite1;
    delete[] fourSite2;

    this->lenX=lenXIn;
    this->lenY=lenYIn;
    this->numSite=lenXIn*lenYIn;

    /*construct bondX and bondY*/
    this->bondX=new int*[this->numSite];
    this->bondY=new int*[this->numSite];
    for(int i=0;i<this->numSite;i++){
        this->bondX[i]=new int[2];
        this->bondY[i]=new int[2];
    }

    int coorTemp[2];
    int codeTemp;
    int bondInd;

    for(int xInd=0;xInd<this->lenX;xInd++){
        for(int yInd=0;yInd<this->lenY;yInd++){
            bondInd=yInd*this->lenX+xInd;
            coorTemp[0]=xInd;
            coorTemp[1]=yInd;
            this->coorToCode(coorTemp,codeTemp);
            this->bondX[bondInd][0]=codeTemp;
            this->bondY[bondInd][0]=codeTemp;

            coorTemp[0]=(xInd+1)%this->lenX;
            coorTemp[1]=yInd;
            this->coorToCode(coorTemp,codeTemp);
            this->bondX[bondInd][1]=codeTemp;

            coorTemp[0]=xInd;
            coorTemp[1]=(yInd+1)%this->lenY;
            this->coorToCode(coorTemp,codeTemp);
            this->bondY[bondInd][1]=codeTemp;
        }
    }

    /*construct fourSite1 and fourSite2*/
    this->fourSite1=new int*[this->numSite/4];
    this->fourSite2=new int*[this->numSite/4];
    for(int i=0;i<this->numSite/4;i++){
        this->fourSite1[i]=new int[4];
        this->fourSite2[i]=new int[4];
    }
    for(int xInd=0;xInd<this->lenX;xInd+=2){
        for(int yInd=0;yInd<this->lenY;yInd+=2){
            int coor0[2]={xInd,yInd};
            int coor1[2]={(xInd+1)%this->lenX,yInd};
            int coor2[2]={(xInd+1)%this->lenX,(yInd+1)%this->lenY};
            int coor3[2]={xInd,(yInd+1)%this->lenY};
            int code0,code1,code2,code3;
            this->coorToCode(coor0,code0);
            this->coorToCode(coor1,code1);
            this->coorToCode(coor2,code2);
            this->coorToCode(coor3,code3);
            this->fourSite1[xInd/2+yInd*this->lenX/4][0]=code0;
            this->fourSite1[xInd/2+yInd*this->lenX/4][1]=code1;
            this->fourSite1[xInd/2+yInd*this->lenX/4][2]=code2;
            this->fourSite1[xInd/2+yInd*this->lenX/4][3]=code3;
        }
    }
    for(int xInd=1;xInd<this->lenX;xInd+=2){
        for(int yInd=1;yInd<this->lenY;yInd+=2){
            int coor0[2]={xInd,yInd};
            int coor1[2]={(xInd+1)%lenX,yInd};
            int coor2[2]={(xInd+1)%lenX,(yInd+1)%this->lenY};
            int coor3[2]={xInd,(yInd+1)%this->lenY};
            int code0,code1,code2,code3;
            this->coorToCode(coor0,code0);
            this->coorToCode(coor1,code1);
            this->coorToCode(coor2,code2);
            this->coorToCode(coor3,code3);
            this->fourSite2[(xInd-1)/2+(yInd-1)*this->lenX/4][0]=code0;
            this->fourSite2[(xInd-1)/2+(yInd-1)*this->lenX/4][1]=code1;
            this->fourSite2[(xInd-1)/2+(yInd-1)*this->lenX/4][2]=code2;
            this->fourSite2[(xInd-1)/2+(yInd-1)*this->lenX/4][3]=code3;
        }
    }
}

void Model::showModelParm(ostream& iout)const{
    iout<<"===============model parameters================"<<endl;
    iout<<"lattice size:            "<<this->lenX<<"x"<<this->lenY<<endl;
    iout<<"hopping strength:        "<<this->strHop<<endl;
    iout<<"hopping ratio:           "<<this->ratHop<<endl;
    iout<<"coupling strength:       "<<this->strCpl<<endl;
    iout<<"dynamic strength:        "<<this->strDyn<<endl;
    iout<<"gradient strength:       "<<this->strGrad<<endl;
    iout<<"squared term strength:   "<<this->strSq<<endl;
    iout<<"quartic term strength:   "<<this->strQt<<endl;
    iout<<"chemical potential:      "<<this->chem<<endl;
    iout<<"chemical potential bias: "<<this->biasChem<<endl;
}



#endif