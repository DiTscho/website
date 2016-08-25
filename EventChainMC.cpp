#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <random>      
#include <functional>  
#include <limits>

const int rng_seed = 42;
// random numbers
std::mt19937 mersenne(rng_seed); // use Mersenne Twister RNG
std::uniform_int_distribution<int> unif(0,3);
std::function<int()> ran03 = std::bind(unif, std::ref(mersenne));

std::uniform_int_distribution<int> unif1(0,1);
std::function<int()> ranBit = std::bind(unif1, std::ref(mersenne));

std::uniform_int_distribution<int> unif2(1,3);
std::function<int()> ran13 = std::bind(unif2, std::ref(mersenne));





using namespace std;


int numbOfdisks = 4;
//int numbOfNextParticles = 1;
double box = 4.0;
double zeroDistanceToGo;

int iterations = 1000000;


vector<double> data;
vector<vector<double> > L;
//vector<double> nextParticle;
double nextParticle[2]={0., 0.};


/*
//periodic distance between two two-dimensional points r1 and r2
double dist(const vector<double>& r1, const vector<double>& r2){
  vector<double> dr(2);
  double drSqd = 0;
  for (int k = 0; k < 2; ++k) {
    dr[k] = r1[k] - r2[k];
    drSqd += dr[k] * dr[k];
}
// if separation in x or y is > box/2 check the closest image
for (int k = 0; k < 2; k++)
    if (sqrt(drSqd) > 0.5 * box)
        dr[k] *= 1 - 1 / sqrt(drSqd);
return sqrt(drSqd);
}
*/
float mod(float a, float N) {return a - N*floor(a/N);}

double dist(const vector<double>& x, const vector<double>& y){
  double d_x =  mod( ( abs(x[0]-y[0]) ), box);
       d_x = min( d_x, box-d_x );
  double d_y = mod( ( abs(x[1]-y[1]) ), box);
       d_y = min(d_y, box-d_y );
  return sqrt(d_x*d_x + d_y*d_y);
}



void flipConf(vector< vector<double> > &v, int &r ){
  vector<vector<double> > Lflip;
  Lflip.resize(v.size());    // numbOfdisks position vectors
  for (int i = 0; i < v.size(); ++i)
       Lflip[i].resize(2, 0);  // each has 2 components

  for (int k = 0; k < v.size(); ++k) {
    if(r==1){
         Lflip[k][0] = box-v[k][1];
         Lflip[k][1] = v[k][0];
    }else{
         Lflip[k][0] = v[k][1];
         Lflip[k][1] = box-v[k][0];
    }
  }
    v = Lflip;
    r = -r;

 }




int main(){

    ofstream outfile("Rij.txt");
    L.resize(numbOfdisks);    // numbOfdisks position vectors
    for (int i = 0; i < numbOfdisks; ++i){L[i].resize(2, 0);}  // each has 2 components
    L = { {1.,1.}, {2.,2.}, {3.,3.}, {3.,1.} };
    int rot = 1;
    std::uniform_real_distribution<double> unif2(0,L.size());
    std::function<int()> choice = std::bind(unif2, std::ref(mersenne));
    double ltilde = 0.9;
for(int iter = 0;iter < iterations;++iter){
        int i = ran03();
        int j = (i+ran13())%4;
        data.push_back( dist( L[i],L[j] ) );
        if (ranBit() < 1) {flipConf(L, rot);}
        zeroDistanceToGo = ltilde;
        int randNum = choice(); 
        nextParticle[0] = L[randNum][0]; 
        nextParticle[1] = L[randNum][1];
    while(zeroDistanceToGo > 0.0){
        L.erase (L.begin()+randNum);
        for (int i = 0; i < L.size(); ++i) { 
			L[i][0] = mod( (L[i][0]-nextParticle[0]), box ); 
			L[i][1] = mod( (L[i][1]-nextParticle[1]), box );
		}
        double nextPosition = std::numeric_limits<double>::infinity(); 
        double a = std::numeric_limits<double>::infinity();
        nextParticle[0] = a; nextParticle[1] = 0.0;
        double currentPosition = 0.0;
        for (int i = 0; i < L.size(); ++i) {
            double xImage[2];
			if (L[i][0]>box/2) {xImage[0] = L[i][0] - box;}
			if (L[i][1]>box/2) {xImage[1] = L[i][1] - box;}
			else{xImage[0]=L[i][0];xImage[1]=L[i][1];}
			if (abs(xImage[1])<1.0) {
				double xDummy = xImage[0] - sqrt( 1.0 -(xImage[1]*xImage[1]) );
				if (xDummy > 0.0 && xDummy < min(zeroDistanceToGo, nextPosition) ) {
					nextPosition = xDummy;
					nextParticle[0] = L[i][0];
					nextParticle[1] = L[i][1];
                }
            }
        }
        double distanceToNextEvent = nextPosition - currentPosition;
        if (zeroDistanceToGo < distanceToNextEvent) {
            currentPosition += zeroDistanceToGo;
            L.push_back( {currentPosition, 0.0} );
        break;
        }else{
            currentPosition += distanceToNextEvent;
            zeroDistanceToGo -= distanceToNextEvent;
            L.push_back( {currentPosition, 0.0} );
            }
    }

}
cout << "----------------------" << endl;
for (int i = 0; i < data.size(); ++i) {
  //cout << std::setprecision(12) << i << "\t" << data[i] << "\n";
  outfile << std::setprecision(12) << data[i] << endl;
}
    return 0;
}
