/*****************************************************
A random class using Mersenne twister algorithm.
Based on mt() in <random>
*****************************************************/
#ifndef _RANDOM_MERSENNE_TWISTER_HPP
#define _RANDOM_MERSENNE_TWISTER_HPP

#include <random>
#include <iostream>
using std::cout;
using std::endl;

class RanMT{
    public:
    RanMT(const int seedIn=0);
    RanMT(const RanMT& obj);
    RanMT& operator =(const RanMT& obj);
    ~RanMT();

    int readSeed()const;
    void setSeed(const int mtSeedIn);
    int ranInt(const int ranMin=0,const int ranMax=1);
    double ranReal(const double ranMin=0,const double ranMax=1);
    int mtSeed;
    std::mt19937 mt;
};

RanMT::RanMT(const int seedIn){
    this->mtSeed=seedIn;
    this->mt.seed(this->mtSeed);
}
RanMT::RanMT(const RanMT& obj){
    this->mtSeed=obj.mtSeed;
    this->mt.seed(this->mtSeed);
}
RanMT& RanMT::operator =(const RanMT& obj){
    if(this==&obj){
        return *this;
    }
    this->mtSeed=obj.mtSeed;
    this->mt.seed(this->mtSeed);
    return *this;
}
RanMT::~RanMT(){

}

int RanMT::readSeed()const{
    return this->mtSeed;
}

void RanMT::setSeed(const int mtSeedIn){
    this->mtSeed=mtSeedIn;
    this->mt.seed(mtSeedIn);
}

int RanMT::ranInt(const int ranMin,const int ranMax){
    std::mt19937::result_type n;
    int range=ranMax-ranMin+1;
    do{
        n=mt();
    }while(n>(mt.max()-(mt.max()-range+1)%range) );
    return n%(ranMax-ranMin+1)+ranMin;
}
double RanMT::ranReal(const double ranMin,const double ranMax){
    double ran01=double(this->mt())/this->mt.max();
    return ranMin+(ranMax-ranMin)*ran01;
}

#endif