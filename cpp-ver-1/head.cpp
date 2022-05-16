#include "head.h"

void Lattice::coorToCode(const int* const coor,int& code){
    code=coor[0]+coor[1]*this->numCellX+coor[2]*this->numCell;
}

void Lattice::codeToCoor(const int code,int* coor){
    coor[2]=code/this->numCell;
    int remain=code-coor[2]*this->numCell;
    coor[1]=remain/this->numCellX;
    coor[0]=remain-coor[1]*this->numCellX;
}

Lattice::Lattice(const int numCellXIn,const int numCellYIn){
    this->numCellX=numCellXIn;
    this->numCellY=numCellYIn;
    this->numCell=this->numCellX*this->numCellY;
    this->numSite=2*this->numCell;

    /*generate the 1st nearest bond lists*/
    this->bondList1X=new int*[this->numCell];
    this->bondList1Y=new int*[this->numCell];
    this->bondList1Z=new int*[this->numCell];
    this->bondList1=new int*[3*this->numCell];
    for(int i=0;i<this->numCell;i++){
        this->bondList1X[i]=new int[2];
        this->bondList1Y[i]=new int[2];
        this->bondList1Z[i]=new int[2];
    }
    for(int i=0;i<3*this->numCell;i++){
        this->bondList1[i]=new int[2];
    }
    for(int xInd=0;xInd<this->numCellX;xInd++){
        for(int yInd=0;yInd<this->numCellY;yInd++){
            int temp1[3]={xInd,(yInd+1)%this->numCellY,0};
            this->coorToCode(temp1,this->bondList1X[yInd*this->numCellX+xInd][0]);
            int temp2[3]={xInd,yInd,1};
            this->coorToCode(temp2,this->bondList1X[yInd*this->numCellX+xInd][1]);
        }
    }
    for(int xInd=0;xInd<this->numCellX;xInd++){
        for(int yInd=0;yInd<this->numCellY;yInd++){
            int temp1[3]={(xInd-1+this->numCellX)%this->numCellX,(yInd+1)%this->numCellY,0};
            this->coorToCode(temp1,this->bondList1Y[yInd*this->numCellX+xInd][0]);
            int temp2[3]={xInd,yInd,1};
            this->coorToCode(temp2,this->bondList1Y[yInd*this->numCellX+xInd][1]);
        }
    }
    for(int xInd=0;xInd<this->numCellX;xInd++){
        for(int yInd=0;yInd<this->numCellY;yInd++){
            int temp1[3]={xInd,yInd,0};
            this->coorToCode(temp1,this->bondList1Z[yInd*this->numCellX+xInd][0]);
            int temp2[3]={xInd,yInd,1};
            this->coorToCode(temp2,this->bondList1Z[yInd*this->numCellX+xInd][1]);
        }
    }
    for(int i=0;i<this->numCell;i++){
        for(int j=0;j<2;j++){
            this->bondList1[i][j]=this->bondList1X[i][j];
            this->bondList1[i+this->numCell][j]=this->bondList1Y[i][j];
            this->bondList1[i+2*this->numCell][j]=this->bondList1Z[i][j];
        }
    }

    /*generate the 2nd nearest bond lists*/
    this->bondList2X=new int*[2*this->numCell];
    this->bondList2Y=new int*[2*this->numCell];
    this->bondList2Z=new int*[2*this->numCell];
    this->bondList2=new int*[6*this->numCell];
    for(int i=0;i<2*this->numCell;i++){
        this->bondList2X[i]=new int[2];
        this->bondList2Y[i]=new int[2];
        this->bondList2Z[i]=new int[2];
    }
    for(int i=0;i<6*this->numCell;i++){
        this->bondList2[i]=new int[2];
    }
    for(int xInd=0;xInd<this->numCellX;xInd++){
        for(int yInd=0;yInd<this->numCellY;yInd++){
            int temp1[3]={xInd,yInd,0};
            this->coorToCode(temp1,this->bondList2X[yInd*this->numCellX+xInd][0]);
            int temp2[3]={(xInd+1)%this->numCellX,(yInd-1+this->numCellY)%this->numCellY,0};
            this->coorToCode(temp2,this->bondList2X[yInd*this->numCellX+xInd][1]);
            temp1[2]=1;
            temp2[2]=1;
            this->coorToCode(temp1,this->bondList2X[yInd*this->numCellX+xInd+this->numCell][0]);
            this->coorToCode(temp2,this->bondList2X[yInd*this->numCellX+xInd+this->numCell][1]);
        }
    }
    for(int xInd=0;xInd<this->numCellX;xInd++){
        for(int yInd=0;yInd<this->numCellY;yInd++){
            int temp1[3]={xInd,yInd,0};
            this->coorToCode(temp1,this->bondList2Y[yInd*this->numCellX+xInd][0]);
            int temp2[3]={xInd,(yInd+1)%this->numCellY,0};
            this->coorToCode(temp2,this->bondList2Y[yInd*this->numCellX+xInd][1]);
            temp1[2]=1;
            temp2[2]=1;
            this->coorToCode(temp1,this->bondList2Y[yInd*this->numCellX+xInd+this->numCell][0]);
            this->coorToCode(temp2,this->bondList2Y[yInd*this->numCellX+xInd+this->numCell][1]);
        }
    }
    for(int xInd=0;xInd<this->numCellX;xInd++){
        for(int yInd=0;yInd<this->numCellY;yInd++){
            int temp1[3]={xInd,yInd,0};
            this->coorToCode(temp1,this->bondList2Z[yInd*this->numCellX+xInd][0]);
            int temp2[3]={(xInd+1)%this->numCellX,yInd,0};
            this->coorToCode(temp2,this->bondList2Z[yInd*this->numCellX+xInd][1]);
            temp1[2]=1;
            temp2[2]=1;
            this->coorToCode(temp1,this->bondList2Z[yInd*this->numCellX+xInd+this->numCell][0]);
            this->coorToCode(temp2,this->bondList2Z[yInd*this->numCellX+xInd+this->numCell][1]);
        }
    }
    for(int i=0;i<2*this->numCell;i++){
        for(int j=0;j<2;j++){
            this->bondList2[i][j]=this->bondList2X[i][j];
            this->bondList2[i+2*this->numCell][j]=bondList2Y[i][j];
            this->bondList2[i+4*this->numCell][j]=bondList2Z[i][j];
        }
    }

    /*generate the 3rd nearest bond lists*/
    this->bondList3X=new int*[this->numCell];
    this->bondList3Y=new int*[this->numCell];
    this->bondList3Z=new int*[this->numCell];
    this->bondList3=new int*[3*this->numCell];
    for(int i=0;i<this->numCell;i++){
        this->bondList3X[i]=new int[2];
        this->bondList3Y[i]=new int[2];
        this->bondList3Z[i]=new int[2];
    }
    for(int i=0;i<3*this->numCell;i++){
        this->bondList3[i]=new int[2];
    }
    for(int xInd=0;xInd<this->numCellX;xInd++){
        for(int yInd=0;yInd<this->numCellY;yInd++){
            int temp1[3]={(xInd-1+this->numCellX)%this->numCellX,yInd,0};
            this->coorToCode(temp1,this->bondList3X[yInd*this->numCellX+xInd][0]);
            int temp2[3]={xInd,yInd,1};
            this->coorToCode(temp2,this->bondList3X[yInd*this->numCellX+xInd][1]);
        }
    }
    for(int xInd=0;xInd<this->numCellX;xInd++){
        for(int yInd=0;yInd<this->numCellY;yInd++){
            int temp1[3]={(xInd+1)%this->numCellX,yInd,0};
            this->coorToCode(temp1,this->bondList3Y[yInd*this->numCellX+xInd][0]);
            int temp2[3]={xInd,yInd,1};
            this->coorToCode(temp2,this->bondList3Y[yInd*this->numCellX+xInd][1]);
        }
    }
    for(int xInd=0;xInd<this->numCellX;xInd++){
        for(int yInd=0;yInd<this->numCellY;yInd++){
            int temp1[3]={(xInd-1+this->numCellX)%this->numCellX,(yInd+2)%this->numCellY,0};
            this->coorToCode(temp1,this->bondList3Z[yInd*this->numCellX+xInd][0]);
            int temp2[3]={xInd,yInd,1};
            this->coorToCode(temp2,this->bondList3Z[yInd*this->numCellX+xInd][1]);
        }
    }
    for(int i=0;i<this->numCell;i++){
        for(int j=0;j<2;j++){
            this->bondList3[i][j]=this->bondList3X[i][j];
            this->bondList3[i+this->numCell][j]=this->bondList3Y[i][j];
            this->bondList3[i+2*this->numCell][j]=this->bondList3Z[i][j];
        }
    }

    /*generate a list: set site code as the index, store 1st nearest neighbors in x y z order*/
    this->siteConnList1=new int*[2*this->numCell];
    for(int i=0;i<2*this->numCell;i++){
        this->siteConnList1[i]=new int[3];
    }
    for(int xInd=0;xInd<this->numCellX;xInd++){
        for(int yInd=0;yInd<this->numCellY;yInd++){
            int code,codeX,codeY,codeZ;
            int coor[3]={xInd,yInd,0};
            int coorX[3]={xInd,(yInd-1+this->numCellY)%this->numCellY,1};
            int coorY[3]={(xInd+1)%this->numCellX,(yInd-1+this->numCellY)%this->numCellY,1};
            int coorZ[3]={xInd,yInd,1};
            this->coorToCode(coor,code);
            this->coorToCode(coorX,codeX);
            this->coorToCode(coorY,codeY);
            this->coorToCode(coorZ,codeZ);
            this->siteConnList1[code][0]=codeX;
            this->siteConnList1[code][1]=codeY;
            this->siteConnList1[code][2]=codeZ;
        }
    }
    for(int xInd=0;xInd<this->numCellX;xInd++){
        for(int yInd=0;yInd<this->numCellY;yInd++){
            int code,codeX,codeY,codeZ;
            int coor[3]={xInd,yInd,1};
            int coorX[3]={xInd,(yInd+1)%this->numCellY,0};
            int coorY[3]={(xInd-1+this->numCellX)%this->numCellX,(yInd+1)%this->numCellY,0};
            int coorZ[3]={xInd,yInd,0};
            this->coorToCode(coor,code);
            this->coorToCode(coorX,codeX);
            this->coorToCode(coorY,codeY);
            this->coorToCode(coorZ,codeZ);
            this->siteConnList1[code][0]=codeX;
            this->siteConnList1[code][1]=codeY;
            this->siteConnList1[code][2]=codeZ;
        }
    }

    /*generate a list: set site code as the index, store 2nd nearest neighbors in xx yy zz order*/
    this->siteConnList2=new int*[2*this->numCell];
    for(int i=0;i<2*this->numCell;i++){
        this->siteConnList2[i]=new int[6];
    }
    for(int xInd=0;xInd<this->numCellX;xInd++){
        for(int yInd=0;yInd<this->numCellY;yInd++){
            int code,codeX1,codeX2,codeY1,codeY2,codeZ1,codeZ2;
            int coor[3]={xInd,yInd,0};
            int coorX1[3]={(xInd+1)%this->numCellX,(yInd-1+this->numCellY)%this->numCellY,0};
            int coorX2[3]={(xInd-1+this->numCellX)%this->numCellX,(yInd+1)%this->numCellY,0};
            int coorY1[3]={xInd,(yInd-1+this->numCellY)%this->numCellY,0};
            int coorY2[3]={xInd,(yInd+1)%this->numCellY,0};
            int coorZ1[3]={(xInd-1+this->numCellX)%this->numCellX,yInd,0};
            int coorZ2[3]={(xInd+1)%this->numCellX,yInd,0};
            this->coorToCode(coor,code);
            this->coorToCode(coorX1,codeX1);
            this->coorToCode(coorX2,codeX2);
            this->coorToCode(coorY1,codeY1);
            this->coorToCode(coorY2,codeY2);
            this->coorToCode(coorZ1,codeZ1);
            this->coorToCode(coorZ2,codeZ2);
            int temp[6]={codeX1,codeX2,codeY1,codeY2,codeZ1,codeZ2};
            for(int i=0;i<6;i++){
                this->siteConnList2[code][i]=temp[i];
            }
        }
    }
    for(int xInd=0;xInd<this->numCellX;xInd++){
        for(int yInd=0;yInd<this->numCellY;yInd++){
            int code,codeX1,codeX2,codeY1,codeY2,codeZ1,codeZ2;
            int coor[3]={xInd,yInd,1};
            int coorX1[3]={(xInd+1)%this->numCellX,(yInd-1+this->numCellY)%this->numCellY,1};
            int coorX2[3]={(xInd-1+this->numCellX)%this->numCellX,(yInd+1)%this->numCellY,1};
            int coorY1[3]={xInd,(yInd-1+this->numCellY)%this->numCellY,1};
            int coorY2[3]={xInd,(yInd+1)%this->numCellY,1};
            int coorZ1[3]={(xInd-1+this->numCellX)%this->numCellX,yInd,1};
            int coorZ2[3]={(xInd+1)%this->numCellX,yInd,1};
            this->coorToCode(coor,code);
            this->coorToCode(coorX1,codeX1);
            this->coorToCode(coorX2,codeX2);
            this->coorToCode(coorY1,codeY1);
            this->coorToCode(coorY2,codeY2);
            this->coorToCode(coorZ1,codeZ1);
            this->coorToCode(coorZ2,codeZ2);
            int temp[6]={codeX1,codeX2,codeY1,codeY2,codeZ1,codeZ2};
            for(int i=0;i<6;i++){
                this->siteConnList2[code][i]=temp[i];
            }
        }
    }

    /*generate a list: set site code as the index, store 3rd nearest neighbors in x y z order*/
    this->siteConnList3=new int*[2*this->numCell];
    for(int i=0;i<2*this->numCell;i++){
        this->siteConnList3[i]=new int[3];
    }
    for(int xInd=0;xInd<this->numCellX;xInd++){
        for(int yInd=0;yInd<this->numCellY;yInd++){
            int code,codeX,codeY,codeZ;
            int coor[3]={xInd,yInd,0};
            int coorX[3]={(xInd+1)%this->numCellX,yInd,1};
            int coorY[3]={(xInd-1+this->numCellX)%this->numCellX,yInd,1};
            int coorZ[3]={(xInd+1)%this->numCellX,(yInd-2+2*this->numCellY)%this->numCellY,1};
            this->coorToCode(coor,code);
            this->coorToCode(coorX,codeX);
            this->coorToCode(coorY,codeY);
            this->coorToCode(coorZ,codeZ);
            this->siteConnList3[code][0]=codeX;
            this->siteConnList3[code][1]=codeY;
            this->siteConnList3[code][2]=codeZ;
        }
    }
    for(int xInd=0;xInd<this->numCellX;xInd++){
        for(int yInd=0;yInd<this->numCellY;yInd++){
            int code,codeX,codeY,codeZ;
            int coor[3]={xInd,yInd,1};
            int coorX[3]={(xInd-1+this->numCellX)%this->numCellX,yInd,0};
            int coorY[3]={(xInd+1)%this->numCellX,yInd,0};
            int coorZ[3]={(xInd-1+this->numCellX)%this->numCellX,(yInd+2)%this->numCellY,0};
            this->coorToCode(coor,code);
            this->coorToCode(coorX,codeX);
            this->coorToCode(coorY,codeY);
            this->coorToCode(coorZ,codeZ);
            this->siteConnList3[code][0]=codeX;
            this->siteConnList3[code][1]=codeY;
            this->siteConnList3[code][2]=codeZ;
        }
    }
}

Lattice::~Lattice(){
    for(int i=0;i<this->numCell;i++){
        delete[] this->bondList1X[i];
        delete[] this->bondList1Y[i];
        delete[] this->bondList1Z[i];
    }
    for(int i=0;i<3*this->numCell;i++){
        delete[] this->bondList1[i];
    }
    delete[] this->bondList1X;
    delete[] this->bondList1Y;
    delete[] this->bondList1Z;
    delete[] this->bondList1;

    for(int i=0;i<2*this->numCell;i++){
        delete[] this->bondList2X[i];
        delete[] this->bondList2Y[i];
        delete[] this->bondList2Z[i];
    }
    for(int i=0;i<6*this->numCell;i++){
        delete[] this->bondList2[i];
    }
    delete[] this->bondList2X;
    delete[] this->bondList2Y;
    delete[] this->bondList2Z;
    delete[] this->bondList2;

    for(int i=0;i<this->numCell;i++){
        delete[] this->bondList3X[i];
        delete[] this->bondList3Y[i];
        delete[] this->bondList3Z[i];
    }
    for(int i=0;i<3*this->numCell;i++){
        delete[] this->bondList3[i];
    }
    delete[] this->bondList3X;
    delete[] this->bondList3Y;
    delete[] this->bondList3Z;
    delete[] this->bondList3;

    for(int i=0;i<2*this->numCell;i++){
        delete[] this->siteConnList1[i];
        delete[] this->siteConnList2[i];
        delete[] this->siteConnList3[i];
    }
    delete this->siteConnList1;
    delete this->siteConnList2;
    delete this->siteConnList3;
}

Lattice::Lattice(const Lattice& obj){
    this->numCellX=obj.numCellX;
    this->numCellY=obj.numCellY;
    this->numCell=obj.numCell;
    this->numSite=obj.numSite;

    /*copy the 1st nearest bond lists*/
    this->bondList1X=new int*[this->numCell];
    this->bondList1Y=new int*[this->numCell];
    this->bondList1Z=new int*[this->numCell];
    this->bondList1=new int*[3*this->numCell];
    for(int i=0;i<this->numCell;i++){
        this->bondList1X[i]=new int[2];
        this->bondList1Y[i]=new int[2];
        this->bondList1Z[i]=new int[2];
        for(int j=0;j<2;j++){
            this->bondList1X[i][j]=obj.bondList1X[i][j];
            this->bondList1Y[i][j]=obj.bondList1Y[i][j];
            this->bondList1Z[i][j]=obj.bondList1Z[i][j];
        }
    }
    for(int i=0;i<3*this->numCell;i++){
        this->bondList1[i]=new int[2];
        for(int j=0;j<2;j++){
            this->bondList1[i][j]=obj.bondList1[i][j];
        }
    }

    /*copy the 2nd nearest bond lists*/
    this->bondList2X=new int*[2*this->numCell];
    this->bondList2Y=new int*[2*this->numCell];
    this->bondList2Z=new int*[2*this->numCell];
    this->bondList2=new int*[6*this->numCell];
    for(int i=0;i<2*this->numCell;i++){
        this->bondList2X[i]=new int[2];
        this->bondList2Y[i]=new int[2];
        this->bondList2Z[i]=new int[2];
        for(int j=0;j<2;j++){
            this->bondList2X[i][j]=obj.bondList2X[i][j];
            this->bondList2Y[i][j]=obj.bondList2Y[i][j];
            this->bondList2Z[i][j]=obj.bondList2Z[i][j];
        }
    }
    for(int i=0;i<6*this->numCell;i++){
        this->bondList2[i]=new int[2];
        for(int j=0;j<2;j++){
            this->bondList2[i][j]=obj.bondList2[i][j];
        }
    }

    /*copy the 3rd nearest bond lists*/
    this->bondList3X=new int*[this->numCell];
    this->bondList3Y=new int*[this->numCell];
    this->bondList3Z=new int*[this->numCell];
    this->bondList3=new int*[3*this->numCell];
    for(int i=0;i<this->numCell;i++){
        this->bondList3X[i]=new int[2];
        this->bondList3Y[i]=new int[2];
        this->bondList3Z[i]=new int[2];
        for(int j=0;j<2;j++){
            this->bondList3X[i][j]=obj.bondList3X[i][j];
            this->bondList3Y[i][j]=obj.bondList3Y[i][j];
            this->bondList3Z[i][j]=obj.bondList3Z[i][j];
        }
    }
    for(int i=0;i<3*this->numCell;i++){
        this->bondList3[i]=new int[2];
        for(int j=0;j<2;j++){
            this->bondList3[i][j]=obj.bondList3[i][j];
        }
    }

    /*copy the 1st nearest connected sites list*/
    this->siteConnList1=new int*[2*this->numCell];
    for(int i=0;i<2*this->numCell;i++){
        this->siteConnList1[i]=new int[3];
        for(int j=0;j<3;j++){
            this->siteConnList1[i][j]=obj.siteConnList1[i][j];
        }
    }

    /*copy the 2nd nearest connected sites list*/
    this->siteConnList2=new int*[2*this->numCell];
    for(int i=0;i<2*this->numCell;i++){
        this->siteConnList2[i]=new int[6];
        for(int j=0;j<6;j++){
            this->siteConnList2[i][j]=obj.siteConnList2[i][j];
        }
    }

    /*copy the 3rd nearest connected sites list*/
    this->siteConnList3=new int*[2*this->numCell];
    for(int i=0;i<2*this->numCell;i++){
        this->siteConnList3[i]=new int[3];
        for(int j=0;j<3;j++){
            this->siteConnList3[i][j]=obj.siteConnList3[i][j];
        }
    }
}

Lattice& Lattice::operator=(const Lattice& obj){
    if(this==&obj){
        return *this;
    }

    for(int i=0;i<this->numCell;i++){
        delete[] this->bondList1X[i];
        delete[] this->bondList1Y[i];
        delete[] this->bondList1Z[i];
    }
    for(int i=0;i<3*this->numCell;i++){
        delete[] this->bondList1[i];
    }
    delete[] this->bondList1X;
    delete[] this->bondList1Y;
    delete[] this->bondList1Z;
    delete[] this->bondList1;

    for(int i=0;i<2*this->numCell;i++){
        delete[] this->bondList2X[i];
        delete[] this->bondList2Y[i];
        delete[] this->bondList2Z[i];
    }
    for(int i=0;i<6*this->numCell;i++){
        delete[] this->bondList2[i];
    }
    delete[] this->bondList2X;
    delete[] this->bondList2Y;
    delete[] this->bondList2Z;
    delete[] this->bondList2;

    for(int i=0;i<this->numCell;i++){
        delete[] this->bondList3X[i];
        delete[] this->bondList3Y[i];
        delete[] this->bondList3Z[i];
    }
    for(int i=0;i<3*this->numCell;i++){
        delete[] this->bondList3[i];
    }
    delete[] this->bondList3X;
    delete[] this->bondList3Y;
    delete[] this->bondList3Z;
    delete[] this->bondList3;

    for(int i=0;i<2*this->numCell;i++){
        delete[] this->siteConnList1[i];
        delete[] this->siteConnList2[i];
        delete[] this->siteConnList3[i];
    }
    delete[] this->siteConnList1;
    delete[] this->siteConnList2;
    delete[] this->siteConnList3;


    this->numCellX=obj.numCellX;
    this->numCellY=obj.numCellY;
    this->numCell=obj.numCell;
    this->numSite=obj.numSite;

    /*copy the 1st nearest bond lists*/
    this->bondList1X=new int*[this->numCell];
    this->bondList1Y=new int*[this->numCell];
    this->bondList1Z=new int*[this->numCell];
    this->bondList1=new int*[3*this->numCell];
    for(int i=0;i<this->numCell;i++){
        this->bondList1X[i]=new int[2];
        this->bondList1Y[i]=new int[2];
        this->bondList1Z[i]=new int[2];
        for(int j=0;j<2;j++){
            this->bondList1X[i][j]=obj.bondList1X[i][j];
            this->bondList1Y[i][j]=obj.bondList1Y[i][j];
            this->bondList1Z[i][j]=obj.bondList1Z[i][j];
        }
    }
    for(int i=0;i<3*this->numCell;i++){
        this->bondList1[i]=new int[2];
        for(int j=0;j<2;j++){
            this->bondList1[i][j]=obj.bondList1[i][j];
        }
    }

    /*copy the 2nd nearest bond lists*/
    this->bondList2X=new int*[2*this->numCell];
    this->bondList2Y=new int*[2*this->numCell];
    this->bondList2Z=new int*[2*this->numCell];
    this->bondList2=new int*[6*this->numCell];
    for(int i=0;i<2*this->numCell;i++){
        this->bondList2X[i]=new int[2];
        this->bondList2Y[i]=new int[2];
        this->bondList2Z[i]=new int[2];
        for(int j=0;j<2;j++){
            this->bondList2X[i][j]=obj.bondList2X[i][j];
            this->bondList2Y[i][j]=obj.bondList2Y[i][j];
            this->bondList2Z[i][j]=obj.bondList2Z[i][j];
        }
    }
    for(int i=0;i<6*this->numCell;i++){
        this->bondList2[i]=new int[2];
        for(int j=0;j<2;j++){
            this->bondList2[i][j]=obj.bondList2[i][j];
        }
    }

    /*copy the 3rd nearest bond lists*/
    this->bondList3X=new int*[this->numCell];
    this->bondList3Y=new int*[this->numCell];
    this->bondList3Z=new int*[this->numCell];
    this->bondList3=new int*[3*this->numCell];
    for(int i=0;i<this->numCell;i++){
        this->bondList3X[i]=new int[2];
        this->bondList3Y[i]=new int[2];
        this->bondList3Z[i]=new int[2];
        for(int j=0;j<2;j++){
            this->bondList3X[i][j]=obj.bondList3X[i][j];
            this->bondList3Y[i][j]=obj.bondList3Y[i][j];
            this->bondList3Z[i][j]=obj.bondList3Z[i][j];
        }
    }
    for(int i=0;i<3*this->numCell;i++){
        this->bondList3[i]=new int[2];
        for(int j=0;j<2;j++){
            this->bondList3[i][j]=obj.bondList3[i][j];
        }
    }

    /*copy the 1st nearest connected sites list*/
    this->siteConnList1=new int*[2*this->numCell];
    for(int i=0;i<2*this->numCell;i++){
        this->siteConnList1[i]=new int[3];
        for(int j=0;j<3;j++){
            this->siteConnList1[i][j]=obj.siteConnList1[i][j];
        }
    }

    /*copy the 2nd nearest connected sites list*/
    this->siteConnList2=new int*[2*this->numCell];
    for(int i=0;i<2*this->numCell;i++){
        this->siteConnList2[i]=new int[6];
        for(int j=0;j<6;j++){
            this->siteConnList2[i][j]=obj.siteConnList2[i][j];
        }
    }

    /*copy the 3rd nearest connected sites list*/
    this->siteConnList3=new int*[2*this->numCell];
    for(int i=0;i<2*this->numCell;i++){
        this->siteConnList3[i]=new int[3];
        for(int j=0;j<3;j++){
            this->siteConnList3[i][j]=obj.siteConnList3[i][j];
        }
    }

    return *this;
}

RanMT::RanMT(const unsigned int mtSeed){
    this->curSeed=mtSeed;
    this->mt.seed(mtSeed);
    this->mtMax=this->mt.max();
}
RanMT::~RanMT(){

}
unsigned int RanMT::readSeed()const{
    return this->curSeed;
}
void RanMT::setSeed(const unsigned int mtSeed){
    this->curSeed=mtSeed;
    this->mt.seed(mtSeed);
}
int RanMT::ranInt(const int ranMin,const int ranMax){
    return this->mt()%(ranMax-ranMin+1);
}
double RanMT::ranReal(const double ranMin,const double ranMax){
    double ran01=double(this->mt())/this->mtMax;
    return ranMin+(ranMax-ranMin)*ran01;
}

SysMC::SysMC(){
    this->invTem=0;
    this->sprAng=PI;
    this->ranMT.setSeed(0);

    this->prmHsb1=0;
    this->prmHsb2=0;
    this->prmHsb3=0;
    this->prmKtv=0;
    this->prmGmm=0;
    this->prmGmmAsym=0;
    this->prmAnIso=0;

    this->conf=new double*[2*lat.numCell];
    for(int i=0;i<2*lat.numCell;i++){
        this->conf[i]=new double[2];
        for(int j=0;j<2;j++){
            this->conf[i][j]=0;
        }
    }
}

SysMC::~SysMC(){
    for(int i=0;i<2*lat.numCell;i++){
        delete[] this->conf[i];
    }
    delete[] this->conf;
}

SysMC::SysMC(const SysMC& obj){
    this->invTem=obj.invTem;
    this->sprAng=obj.sprAng;
    this->ranMT.setSeed(obj.ranMT.readSeed());

    this->prmHsb1=obj.prmHsb1;
    this->prmHsb2=obj.prmHsb2;
    this->prmHsb3=obj.prmHsb3;
    this->prmKtv=obj.prmKtv;
    this->prmGmm=obj.prmGmm;
    this->prmGmmAsym=obj.prmGmmAsym;

    this->conf=new double*[2*lat.numCell];
    for(int i=0;i<2*lat.numCell;i++){
        this->conf[i]=new double[2];
        for(int j=0;j<2;j++){
            this->conf[i][j]=obj.conf[i][j];
        }
    }
}

SysMC& SysMC::operator=(const SysMC& obj){
    if(this==&obj){
        return *this;
    }
    this->invTem=obj.invTem;
    this->sprAng=obj.sprAng;
    this->ranMT.setSeed(obj.ranMT.readSeed());

    this->prmHsb1=obj.prmHsb1;
    this->prmHsb2=obj.prmHsb2;
    this->prmHsb3=obj.prmHsb3;
    this->prmKtv=obj.prmKtv;
    this->prmGmm=obj.prmGmm;
    this->prmGmmAsym=obj.prmGmmAsym;

    for(int i=0;i<2*lat.numCell;i++){
        for(int j=0;j<2;j++){
            this->conf[i][j]=obj.conf[i][j];
        }
    }
    return *this;
}

void SysMC::showConf(ostream& dest)const{
    for(int i=0;i<2*lat.numCell;i++){
        for(int j=0;j<2;j++){
            dest<<this->conf[i][j]<<" ";
        }
        dest<<endl;
    }
}

void SysMC::showConfCar(ostream& dest)const{
    double** confCar=new double*[2*lat.numCell];
    for(int i=0;i<2*lat.numCell;i++){
        confCar[i]=new double[3];
        sphToCar(this->conf[i],confCar[i]);
        for(int j=0;j<3;j++){
            dest<<confCar[i][j]<<" ";
        }
        dest<<endl;
    }

    for(int i=0;i<2*lat.numCell;i++){
        delete[] confCar[i];
    }
    delete[] confCar;
}

void SysMC::showPrmMC(ostream& dest)const{
    dest<<"===============Monte Carlo parameters=================="<<endl;
    dest<<"inverse temperature: "<<this->invTem<<endl;
    dest<<"update spread angle: "<<this->sprAng<<endl;
    dest<<"random seed:         "<<this->ranMT.readSeed()<<endl;
    dest<<"======================================================="<<endl;
 }

void SysMC::showPrmHmlt(ostream& dest)const{
    dest<<"===============Hamiltonian parameters=================="<<endl;
    dest<<"Heisenberg (1st nearest neighbor): "<<this->prmHsb1<<endl;
    dest<<"Heisenberg (2nd nearest neighbor): "<<this->prmHsb2<<endl;
    dest<<"Heisenberg (3rd nearest neighbor): "<<this->prmHsb3<<endl;
    dest<<"Kitaev:                            "<<this->prmKtv<<endl;
    dest<<"Symmetric off-diagonal (Gamma):    "<<this->prmGmm<<endl;
    dest<<"Asymmetric off-diagonal (Gamma):   "<<this->prmGmmAsym<<endl;
    dest<<"Anisotropic (easy-plane):          "<<this->prmAnIso<<endl;
    dest<<"======================================================="<<endl;
}

void SysMC::setPrmMC(const double invTemIn,const double sprAngIn,const unsigned mtSeedIn){
    this->invTem=invTemIn;
    this->sprAng=sprAngIn;
    this->ranMT.setSeed(mtSeedIn);
}

void SysMC::setPrmHmlt(const vector<double>& prmHmlt,ostream& dest){
    if(prmHmlt.size()!=7){
        dest<<"Number of input parameters not compatible!"<<endl;
        dest<<"Hamiltonian parameters setting failed."<<endl;
        dest<<"Procedure terminated."<<endl;
        return;
    }

    this->prmHsb1=prmHmlt[0];
    this->prmHsb2=prmHmlt[1];
    this->prmHsb3=prmHmlt[2];
    this->prmKtv=prmHmlt[3];
    this->prmGmm=prmHmlt[4];
    this->prmGmmAsym=prmHmlt[5];
    this->prmAnIso=prmHmlt[6];
}

int SysMC::ranSite(){
    return this->ranMT.ranInt(0,2*lat.numCell-1);
}

void SysMC::ranStateSph(double* state){
    state[0]=acos(-2*this->ranMT.ranReal(-0.5,0.5));
    state[1]=this->ranMT.ranReal(0,2*PI);
}

void SysMC::ranStateCone(double* state,const double polCen,const double aziCen,const double sprAng){
    double polTemp=acos(-2*this->ranMT.ranReal(-0.5,-0.5*cos(this->sprAng)));
    double aziTemp=this->ranMT.ranReal(0,2*PI);

    double cosPolCen=cos(polCen);
    double sinPolCen=sin(polCen);
    double cosAziCen=cos(aziCen);
    double sinAziCen=sin(aziCen);

    double cosPolTemp=cos(polTemp);
    double sinPolTemp=sin(polTemp);
    double cosAziTemp=cos(aziTemp);
    double sinAziTemp=sin(aziTemp);


    double a=cosPolCen*cosAziTemp*cosAziCen*sinPolTemp+cosPolTemp*cosAziCen*sinPolCen-sinPolTemp*sinAziTemp*sinAziCen;
    double b=cosAziCen*sinPolTemp*sinAziTemp+cosPolCen*cosAziTemp*sinPolTemp*sinAziCen+cosPolTemp*sinPolCen*sinAziCen;
    double c=cosPolTemp*cosPolCen-cosAziTemp*sinPolTemp*sinPolCen;

    state[0]=acos(c);
    if(b>0){
        state[1]=acos(a/sqrt(1-c*c));
    }else{
        state[1]=2*PI-acos(a/sqrt(1-c*c));
    }
}

void SysMC::randomizeConf(){
    for(int i=0;i<2*lat.numCell;i++){
        double temp[2];
        this->ranStateSph(conf[i]);
    }
}

void SysMC::setConf(double** confIn){
    for(int i=0;i<2*lat.numCell;i++){
        for(int j=0;j<2;j++){
            this->conf[i][j]=confIn[i][j];
        }
    }
}

double SysMC::eneTotHsb1(double** confCar){
    double sum=0;
    for(int bondInd=0;bondInd<3*lat.numCell;bondInd++){
        int siteI=lat.bondList1[bondInd][0];
        int siteJ=lat.bondList1[bondInd][1];
        for(int i=0;i<3;i++){
            sum+=confCar[siteI][i]*confCar[siteJ][i];
        }
    }
    return this->prmHsb1*sum;
}

double SysMC::eneTotHsb2(double** confCar){
    double sum=0;
    for(int bondInd=0;bondInd<6*lat.numCell;bondInd++){
        int siteI=lat.bondList2[bondInd][0];
        int siteJ=lat.bondList2[bondInd][1];
        for(int i=0;i<3;i++){
            sum+=confCar[siteI][i]*confCar[siteJ][i];
        }
    }
    return this->prmHsb2*sum;
}

double SysMC::eneTotHsb3(double** confCar){
    double sum=0;
    for(int bondInd=0;bondInd<3*lat.numCell;bondInd++){
        int siteI=lat.bondList3[bondInd][0];
        int siteJ=lat.bondList3[bondInd][1];
        for(int i=0;i<3;i++){
            sum+=confCar[siteI][i]*confCar[siteJ][i];
        }
    }
    return this->prmHsb3*sum;
}

double SysMC::eneTotKtv(double** confCar){
    double sum=0;
    for(int bondInd=0;bondInd<lat.numCell;bondInd++){
        int siteI=lat.bondList1X[bondInd][0];
        int siteJ=lat.bondList1X[bondInd][1];
        sum+=confCar[siteI][0]*confCar[siteJ][0];
    }
    for(int bondInd=0;bondInd<lat.numCell;bondInd++){
        int siteI=lat.bondList1Y[bondInd][0];
        int siteJ=lat.bondList1Y[bondInd][1];
        sum+=confCar[siteI][1]*confCar[siteJ][1];
    }
    for(int bondInd=0;bondInd<lat.numCell;bondInd++){
        int siteI=lat.bondList1Z[bondInd][0];
        int siteJ=lat.bondList1Z[bondInd][1];
        sum+=confCar[siteI][2]*confCar[siteJ][2];
    }
    return this->prmKtv*sum;
}

double SysMC::eneTotGmm(double** confCar){
    double sum=0;
    for(int bondInd=0;bondInd<lat.numCell;bondInd++){
        int siteI=lat.bondList1X[bondInd][0];
        int siteJ=lat.bondList1X[bondInd][1];
        sum+=confCar[siteI][1]*confCar[siteJ][2];
        sum+=confCar[siteI][2]*confCar[siteJ][1];
    }
    for(int bondInd=0;bondInd<lat.numCell;bondInd++){
        int siteI=lat.bondList1Y[bondInd][0];
        int siteJ=lat.bondList1Y[bondInd][1];
        sum+=confCar[siteI][2]*confCar[siteJ][0];
        sum+=confCar[siteI][0]*confCar[siteJ][2];
    }
    for(int bondInd=0;bondInd<lat.numCell;bondInd++){
        int siteI=lat.bondList1Z[bondInd][0];
        int siteJ=lat.bondList1Z[bondInd][1];
        sum+=confCar[siteI][0]*confCar[siteJ][1];
        sum+=confCar[siteI][1]*confCar[siteJ][0];
    }
    return this->prmGmm*sum;
}

double SysMC::eneTotGmmAsym(double** confCar){
    double sum=0;
    for(int bondInd=0;bondInd<lat.numCell;bondInd++){
        int siteI=lat.bondList1X[bondInd][0];
        int siteJ=lat.bondList1X[bondInd][1];
        sum+=confCar[siteI][1]*confCar[siteJ][2];
        sum-=confCar[siteI][2]*confCar[siteJ][1];
    }
    for(int bondInd=0;bondInd<lat.numCell;bondInd++){
        int siteI=lat.bondList1Y[bondInd][0];
        int siteJ=lat.bondList1Y[bondInd][1];
        sum+=confCar[siteI][2]*confCar[siteJ][0];
        sum-=confCar[siteI][0]*confCar[siteJ][2];
    }
    for(int bondInd=0;bondInd<lat.numCell;bondInd++){
        int siteI=lat.bondList1Z[bondInd][0];
        int siteJ=lat.bondList1Z[bondInd][1];
        sum+=confCar[siteI][0]*confCar[siteJ][1];
        sum-=confCar[siteI][1]*confCar[siteJ][0];
    }
    return this->prmGmmAsym*sum;
}

double SysMC::eneTotAnIso(double** confCar){
    double sum=0;
    for(int siteInd=0;siteInd<2*lat.numCell;siteInd++){
        double temp=0;
        for(int i=0;i<3;i++){
            temp+=confCar[siteInd][i];
        }
        sum+=temp*temp;
    }
    return this->prmAnIso*sum;
}

double SysMC::eneTot(){
    double res=0;
    double** confCar=new double*[2*lat.numCell];
    for(int i=0;i<2*lat.numCell;i++){
        confCar[i]=new double[3];
        sphToCar(this->conf[i],confCar[i]);
    }

    res+=this->eneTotHsb1(confCar);
    res+=this->eneTotHsb2(confCar);
    res+=this->eneTotHsb3(confCar);
    res+=this->eneTotKtv(confCar);
    res+=this->eneTotGmm(confCar);
    res+=this->eneTotGmmAsym(confCar);
    res+=this->eneTotAnIso(confCar);

    for(int i=0;i<2*lat.numCell;i++){
        delete[] confCar[i];
    }
    delete[] confCar;

    return res;
}

double SysMC::eneDiff(const int siteUpd,double* spinUpd){
    int site1X=lat.siteConnList1[siteUpd][0];
    int site1Y=lat.siteConnList1[siteUpd][1];
    int site1Z=lat.siteConnList1[siteUpd][2];

    int site2X1=lat.siteConnList2[siteUpd][0];
    int site2X2=lat.siteConnList2[siteUpd][1];
    int site2Y1=lat.siteConnList2[siteUpd][2];
    int site2Y2=lat.siteConnList2[siteUpd][3];
    int site2Z1=lat.siteConnList2[siteUpd][4];
    int site2Z2=lat.siteConnList2[siteUpd][5];

    int site3X=lat.siteConnList3[siteUpd][0];
    int site3Y=lat.siteConnList3[siteUpd][1];
    int site3Z=lat.siteConnList3[siteUpd][2];

    double spin1X[3];
    double spin1Y[3];
    double spin1Z[3];
    sphToCar(this->conf[site1X],spin1X);
    sphToCar(this->conf[site1Y],spin1Y);
    sphToCar(this->conf[site1Z],spin1Z);

    double spin2X1[3];
    double spin2X2[3];
    double spin2Y1[3];
    double spin2Y2[3];
    double spin2Z1[3];
    double spin2Z2[3];
    sphToCar(this->conf[site2X1],spin2X1);
    sphToCar(this->conf[site2X2],spin2X2);
    sphToCar(this->conf[site2Y1],spin2Y1);
    sphToCar(this->conf[site2Y2],spin2Y2);
    sphToCar(this->conf[site2Z1],spin2Z1);
    sphToCar(this->conf[site2Z2],spin2Z2);

    double spin3X[3];
    double spin3Y[3];
    double spin3Z[3];
    sphToCar(this->conf[site3X],spin3X);
    sphToCar(this->conf[site3Y],spin3Y);
    sphToCar(this->conf[site3Z],spin3Z);

    double spinCurr[3];
    sphToCar(this->conf[siteUpd],spinCurr);
    double spinUpdCar[3];
    sphToCar(spinUpd,spinUpdCar);
    double spinDiff[3];
    for(int i=0;i<3;i++){
        spinDiff[i]=spinUpdCar[i]-spinCurr[i];
    }

    double res=0;

    double temp=0;
    for(int i=0;i<3;i++){
        temp+=spinDiff[i]*spin1X[i];
    }
    for(int i=0;i<3;i++){
        temp+=spinDiff[i]*spin1Y[i];
    }
    for(int i=0;i<3;i++){
        temp+=spinDiff[i]*spin1Z[i];
    }
    res+=this->prmHsb1*temp;

    temp=0;
    for(int i=0;i<3;i++){
        temp+=spinDiff[i]*spin2X1[i];
    }
    for(int i=0;i<3;i++){
        temp+=spinDiff[i]*spin2X2[i];
    }
    for(int i=0;i<3;i++){
        temp+=spinDiff[i]*spin2Y1[i];
    }
    for(int i=0;i<3;i++){
        temp+=spinDiff[i]*spin2Y2[i];
    }
    for(int i=0;i<3;i++){
        temp+=spinDiff[i]*spin2Z1[i];
    }
    for(int i=0;i<3;i++){
        temp+=spinDiff[i]*spin2Z2[i];
    }
    res+=this->prmHsb2*temp;

    temp=0;
    for(int i=0;i<3;i++){
        temp+=spinDiff[i]*spin3X[i];
    }
    for(int i=0;i<3;i++){
        temp+=spinDiff[i]*spin3Y[i];
    }
    for(int i=0;i<3;i++){
        temp+=spinDiff[i]*spin3Z[i];
    }
    res+=this->prmHsb3*temp;

    temp=0;
    temp+=spinDiff[0]*spin1X[0];
    temp+=spinDiff[1]*spin1Y[1];
    temp+=spinDiff[2]*spin1Z[2];
    res+=this->prmKtv*temp;

    temp=0;
    temp+=spinDiff[1]*spin1X[2]+spinDiff[2]*spin1X[1];
    temp+=spinDiff[2]*spin1Y[0]+spinDiff[0]*spin1Y[2];
    temp+=spinDiff[0]*spin1Z[1]+spinDiff[1]*spin1Z[0];
    res+=this->prmGmm*temp;

    temp=0;
    temp+=spinDiff[1]*spin1X[2]-spinDiff[2]*spin1X[1];
    temp+=spinDiff[2]*spin1Y[0]-spinDiff[0]*spin1Y[2];
    temp+=spinDiff[0]*spin1Z[1]-spinDiff[1]*spin1Z[0];
    if(siteUpd<lat.numCell){
        res+=this->prmGmmAsym*temp;
    }else{
        res-=this->prmGmmAsym*temp;
    }

    temp=0;
    for(int i=0;i<3;i++){
        temp+=spinUpdCar[i];
    }
    res+=this->prmAnIso*temp*temp;
    temp=0;
    for(int i=0;i<3;i++){
        temp+=spinCurr[i];
    }
    res-=this->prmAnIso*temp*temp;

    return res;

}

void SysMC::tryUpdSpin(const int siteUpd,double* spinUpd){
    double eneChg=this->eneDiff(siteUpd,spinUpd);
    double accProb=0;
    if(eneChg<0){
        accProb=1;
    }else{
        accProb=exp(-invTem*eneChg);
    }
    if(this->ranMT.ranReal()<accProb){
        this->conf[siteUpd][0]=spinUpd[0];
        this->conf[siteUpd][1]=spinUpd[1];
    }
}

void SysMC::tryUpdSpin(const int siteUpd,double* spinUpd,int& accFlag){
    double eneChg=this->eneDiff(siteUpd,spinUpd);
    double accProb=0;
    if(eneChg<0){
        accProb=1;
    }else{
        accProb=exp(-invTem*eneChg);
    }
    accFlag=0;
    if(this->ranMT.ranReal()<accProb){
        this->conf[siteUpd][0]=spinUpd[0];
        this->conf[siteUpd][1]=spinUpd[1];
        accFlag=1;
    }
}

void SysMC::tryUpdSpin(int& accFlag){
    int siteUpd=this->ranSite();
    double spinUpd[2];
    this->ranStateCone(spinUpd,this->conf[siteUpd][0],this->conf[siteUpd][1],this->sprAng);
    double eneChg=this->eneDiff(siteUpd,spinUpd);
    double accProb=0;
    if(eneChg<0){
        accProb=1;
    }else{
        accProb=exp(-invTem*eneChg);
    }
    accFlag=0;
    if(this->ranMT.ranReal()<accProb){
        this->conf[siteUpd][0]=spinUpd[0];
        this->conf[siteUpd][1]=spinUpd[1];
        accFlag=1;
    }
}

void SysMC::swpRan(){
    for(int trlInd=0;trlInd<2*lat.numCell;trlInd++){
        int siteInd=this->ranSite();
        double spinNew[2];
        this->ranStateCone(spinNew,this->conf[siteInd][0],this->conf[siteInd][1],this->sprAng);
        this->tryUpdSpin(siteInd,spinNew);
    }
}

void SysMC::swpRan(double& accRate){
    int accCnt=0;
    int accFlag=0;
    for(int trlInd=0;trlInd<2*lat.numCell;trlInd++){
        int siteInd=this->ranSite();
        double spinNew[2];
        this->ranStateCone(spinNew,this->conf[siteInd][0],this->conf[siteInd][1],this->sprAng);
        this->tryUpdSpin(siteInd,spinNew,accFlag);
        accCnt+=accFlag;
    }
    accRate=double(accCnt)/2/lat.numCell;
}

void SysMC::swpSqn(){
    for(int siteInd=0;siteInd<2*lat.numCell;siteInd++){
        double spinNew[2];
        this->ranStateCone(spinNew,this->conf[siteInd][0],this->conf[siteInd][1],this->sprAng);
        this->tryUpdSpin(siteInd,spinNew);
    }
}

void SysMC::swpSqn(double& accRate){
    int accCnt=0;
    int accFlag=0;
    for(int siteInd=0;siteInd<2*lat.numCell;siteInd++){
        double spinNew[2];
        this->ranStateCone(spinNew,this->conf[siteInd][0],this->conf[siteInd][1],this->sprAng);
        this->tryUpdSpin(siteInd,spinNew,accFlag);
        accCnt+=accFlag;
    }
    accRate=double(accCnt)/2/lat.numCell;
}

double SysMC::calcAccRateFast(){
    vector<double> arPath(0);
    double arTemp;
    for(int genInd=0;genInd<10000;genInd++){
        this->swpRan(arTemp);
        arPath.push_back(arTemp);
        if((genInd+1)%5==0){
            if(calcErrBin(arPath,5)<0.03){
                break;
            }
        }
    }
    return calcMean(arPath);
}

double SysMC::calcAccRateRan(){
    vector<double> arPath(0);
    double arTemp;
    for(int genInd=0;genInd<10000;genInd++){
        this->swpRan(arTemp);
        arPath.push_back(arTemp);
        if((genInd+1)%5==0&&calcErrBin(arPath,5)<0.005){
            break;
        }
    }
    return calcMean(arPath);
}

double SysMC::calcAccRateSqn(){
    vector<double> arPath(0);
    double arTemp;
    for(int genInd=0;genInd<10000;genInd++){
        this->swpSqn(arTemp);
        arPath.push_back(arTemp);
        if((genInd+1)%5==0&&calcErrBin(arPath,5)<0.005){
            break;
        }
    }
    return calcMean(arPath);
}

void SysMC::adjustWinFast(){
    double accRateCurr=this->calcAccRateFast();
    if( (accRateCurr>0.45)&&(accRateCurr<0.55) ){
        return;
    }

    double winL=0;
    double winR=PI;
    double winMid;
    double accRateCur;
    
    this->sprAng=winL;
    const double accRateL=this->calcAccRateFast();
    this->sprAng=winR;
    const double accRateR=this->calcAccRateFast();
    
    if(accRateR>0.55){
        return;
    }

    for(int genInd=0;genInd<20;genInd++){
        winMid=(winL+winR)/2;
        this->sprAng=winMid;
        accRateCur=this->calcAccRateFast();
        if(accRateCur<=0.45){
            winR=winMid;
        }
        if(accRateCur>=0.55){
            winL=winMid;
        }
        if(accRateCur>0.45&&accRateCur<0.55){
            break;
        }
    }
}

void SysMC::adjustWinRan(ostream& iout){
    double accRateCurr=this->calcAccRateRan();
    if( (accRateCurr>0.45)&&(accRateCurr<0.55) ){
        iout<<"current accept rate is "<<accRateCurr<<endl;
        iout<<"random window adjustment not needed"<<endl;
        return;
    }

    double winL=0;
    double winR=PI;
    double winMid;
    double accRateCur;
    
    this->sprAng=winL;
    const double accRateL=this->calcAccRateRan();
    this->sprAng=winR;
    const double accRateR=this->calcAccRateRan();
    iout<<"accept rate for the left bound : "<<accRateL<<endl;
    iout<<"accept rate for the right bound: "<<accRateR<<endl;
    
    if(accRateR>0.49){
        iout<<"Random window is fixed at "<<this->sprAng<<" with average accept rate "<<accRateR<<"."<<endl;
        return;
    }

    for(int genInd=0;genInd<100;genInd++){
        winMid=(winL+winR)/2;
        this->sprAng=winMid;
        accRateCur=this->calcAccRateRan();
        if(accRateCur<=0.45){
            winR=winMid;
        }
        if(accRateCur>=0.55){
            winL=winMid;
        }
        iout<<"trial window: "<<this->sprAng<<"     accept rate: "<<accRateCur<<endl;
        if(accRateCur>0.45&&accRateCur<0.55){
            break;
        }
    }
    iout<<"Random window is fixed at "<<this->sprAng<<" with average accept rate "<<accRateCur<<"."<<endl;
}

void SysMC::adjustWinSqn(ostream& iout){
    double accRateCurr=this->calcAccRateSqn();
    if( (accRateCurr>0.45)&&(accRateCurr<0.55) ){
        iout<<"current accept rate is "<<accRateCurr<<endl;
        iout<<"random window adjustment not needed"<<endl;
        return;
    }

    double winL=0;
    double winR=PI;
    double winMid;
    double accRateCur;
    
    this->sprAng=winL;
    const double accRateL=this->calcAccRateSqn();
    this->sprAng=winR;
    const double accRateR=this->calcAccRateSqn();
    iout<<"accept rate for the left bound : "<<accRateL<<endl;
    iout<<"accept rate for the right bound: "<<accRateR<<endl;
    
    if(accRateR>0.49){
        iout<<"Random window is fixed at "<<this->sprAng<<" with average accept rate "<<accRateR<<"."<<endl;
        return;
    }

    for(int genInd=0;genInd<100;genInd++){
        winMid=(winL+winR)/2;
        this->sprAng=winMid;
        accRateCur=this->calcAccRateSqn();
        if(accRateCur<=0.45){
            winR=winMid;
        }
        if(accRateCur>=0.55){
            winL=winMid;
        }
        iout<<"trial window: "<<this->sprAng<<"     accept rate: "<<accRateCur<<endl;
        if(accRateCur>0.45&&accRateCur<0.55){
            break;
        }
    }
    iout<<"Random window is fixed at "<<this->sprAng<<" with average accept rate "<<accRateCur<<"."<<endl;
}

void SysMC::calcStrFct(double** res){
    double** confCar=new double*[2*lat.numCell];
    for(int i=0;i<2*lat.numCell;i++){
        confCar[i]=new double[3];
        sphToCar(this->conf[i],confCar[i]);
    }

    /*generate rList*/
    double uniX[2]={1,0};
    double uniY[2]={1.0/2,sqrt(3)/2};
    double uniC[2]={0,1.0/sqrt(3)};
    double** rList=new double*[2*lat.numCell];
    for(int i=0;i<2*lat.numCell;i++){
        rList[i]=new double[2];
    }
    for(int x=0;x<lat.numCellX;x++){
        for(int y=0;y<lat.numCellY;y++){
            int i=y*lat.numCellX+x;
            rList[i][0]=x*uniX[0]+y*uniY[0];
            rList[i][1]=x*uniX[1]+y*uniY[1];
            rList[i+lat.numCell][0]=rList[i][0];
            rList[i+lat.numCell][1]=rList[i][1]+1.0/sqrt(3);
        }
    }

    /*generate kList*/
    double preFacX=4*sqrt(3)*PI/3/lat.numCellX;
    double preFacY=4*sqrt(3)*PI/3/lat.numCellY;
    double uniKX[2]={preFacX*sqrt(3)/2,-preFacX/2};
    double uniKY[2]={0,preFacY};
    double** kList=new double*[lat.numCell];
    for(int i=0;i<lat.numCell;i++){
        kList[i]=new double[2];
    }
    for(int x=0;x<lat.numCellX;x++){
        for(int y=0;y<lat.numCellY;y++){
            int site=y*lat.numCellX+x;
            for(int i=0;i<2;i++){
                kList[site][i]=x*uniKX[i]+y*uniKY[i];
                res[site][i]=kList[site][i];
            }
        }
    }

    /*calulate static structure factor*/
    for(int kInd=0;kInd<lat.numCell;kInd++){
        double sum=0;
        for(int i=0;i<2*lat.numCell;i++){
            for(int j=0;j<2*lat.numCell;j++){
                double term1=0;
                for(int q=0;q<3;q++){
                    term1+=confCar[i][q]*confCar[j][q];
                }
                double term2=0;
                for(int q=0;q<2;q++){
                    term2+=kList[kInd][q]*(rList[i][q]-rList[j][q]);
                }
                sum+=term1*cos(term2);
            }
        }
        res[kInd][2]=sum/2/lat.numCell;
    }

    for(int i=0;i<2*lat.numCell;i++){
        delete[] rList[i];
    }
    delete[] rList;
    for(int i=0;i<lat.numCell;i++){
        delete[] kList[i];
    }
    delete[] kList;
    for(int i=0;i<2*lat.numCell;i++){
        delete[] confCar[i];
    }
    delete[] confCar;
}

void SysMC::calcStrFct(ostream& iout){
    double** strFct=new double*[lat.numCell];
    for(int i=0;i<lat.numCell;i++){
        strFct[i]=new double[3];
    }

    this->calcStrFct(strFct);
    for(int i=0;i<lat.numCell;i++){
        for(int j=0;j<3;j++){
            iout<<strFct[i][j]<<" ";
        }
        iout<<endl;
    }

    for(int i=0;i<lat.numCell;i++){
        delete[] strFct[i];
    }
    delete[] strFct;
}

void SysMC::smltAnnl(vector<double>& invTemPath,const int numSwp,ostream& iout){
    int numInvTem=invTemPath.size();
    for(int invTemInd=0;invTemInd<numInvTem;invTemInd++){
        this->invTem=invTemPath[invTemInd];
        iout<<"step "<<invTemInd<<" with inverse temperature "<<this->invTem<<endl;
        this->adjustWinRan(iout);
        // this->adjustWinSqn(iout);
        for(int swpInd=0;swpInd<numSwp;swpInd++){
            this->swpRan();
            // this->swpSqn();
        }
        iout<<"current energy: "<<this->eneTot()<<endl<<endl;;
    }
}

void SysMC::smltAnnl(const double temIni,const double temEnd,const double sclFct,const int winAdjIntv,const int mntrIntv,ostream& iout){
    clock_t c1=clock();
    clock_t c2;
    
    double invTemIni=1/temIni;
    double invTemEnd=1/temEnd;

    this->invTem=invTemIni;
    double costTot=0;
    double entrTot=0;

    double eneIni=this->eneTot();
    cout<<eneTot()<<endl;

    int numUpdMax=int(pow(10,16));
    for(int updInd=0;updInd<numUpdMax;updInd++){
        // if(updInd<10100&&updInd>9900){
        //     cout<<updInd<<" "<<this->eneTot()<<endl;
        // }
        int siteUpd=this->ranSite();
        double spinUpd[2];
        this->ranStateCone(spinUpd,this->conf[siteUpd][0],this->conf[siteUpd][1],this->sprAng);
        double eneChg=this->eneDiff(siteUpd,spinUpd);
        double accProb=0;
        if(eneChg<0){
            accProb=1;
        }else{
            accProb=exp(-invTem*eneChg);
        }
        if(this->ranMT.ranReal()<accProb){
            this->conf[siteUpd][0]=spinUpd[0];
            this->conf[siteUpd][1]=spinUpd[1];
            costTot+=eneChg;
            if(eneChg>0){
                entrTot-=eneChg*this->invTem;
            }
            if(costTot>=0||entrTot==0){
                this->invTem=invTemIni;
            }else{
                this->invTem=sclFct*entrTot/costTot;
            }
        }

        if( (updInd+1)%winAdjIntv==0 ){
            this->adjustWinFast();
            costTot=this->eneTot()-eneIni;
        }

        if((updInd+1)%mntrIntv==0){
            iout<<updInd+1<<" "<<this->invTem<<" "<<eneIni+costTot<<endl;
        }

        if(updInd>50000&&this->invTem>invTemEnd){
            c2=clock();
            iout<<"annealing procedure finished with "<<updInd+1<<" steps"<<endl;
            iout<<"current inverse temperature: "<<this->invTem<<endl;
            iout<<"current energy:              "<<this->eneTot()<<endl;
            iout<<"time spent on annealing:     "<<double(c2-c1)/CLOCKS_PER_SEC<<" seconds"<<endl;
            return;
        }
    }
}