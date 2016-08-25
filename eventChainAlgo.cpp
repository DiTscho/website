/* no description here, just an implementation of python
   code by Krauth to go sure that further implementations
   will be correct */

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

const int rng_seed = 42;
const int N = 4;
// random numbers
std::mt19937 mersenne(rng_seed); // use Mersenne Twister RNG
std::uniform_int_distribution<int> unif(0, N - 1);
std::function<int()> ran0_Nmin1 = std::bind(unif, std::ref(mersenne));

std::uniform_int_distribution<int> unif1(0, 1);
std::function<int()> ranBit = std::bind(unif1, std::ref(mersenne));

std::uniform_int_distribution<int> unif2(1, N - 1);
std::function<int()> ran1_Nmin1 = std::bind(unif2, std::ref(mersenne));

std::uniform_real_distribution<double> unif4(0, 1);
std::function<double()> ranReal = std::bind(unif4, std::ref(mersenne));

using namespace std;

int iterations = 1000000;
double box = 1.0;
double d = 0.25;
double zeroDistanceToGo;

vector<double> data;
vector<vector<double>> L;
vector<double> lambdaInterim;
vector<double> lambda;
double counts = 0.;
double lambda0 = 0.;
// vector<double> nextParticle;
double nextParticle[2] = {0., 0.};

float mod(float a, float b) { return a - b * floor(a / b); }

double dist(const vector<double> &x, const vector<double> &y) {
  double d_x = mod((abs(x[0] - y[0])), box);
  d_x = min(d_x, box - d_x);
  double d_y = mod((abs(x[1] - y[1])), box);
  d_y = min(d_y, box - d_y);
  return sqrt(d_x * d_x + d_y * d_y);
}

void flipConf(vector<vector<double>> &v, int &r) {
  vector<vector<double>> Lflip;
  Lflip.resize(v.size()); // N position vectors
  for (int i = 0; i < v.size(); ++i)
    Lflip[i].resize(2, 0); // each has 2 components
                           //  Lflip = v;

  for (int k = 0; k < v.size(); ++k) {
    if (r == 1) {
      Lflip[k][0] = box - v[k][1];
      Lflip[k][1] = v[k][0];
    } else {
      Lflip[k][0] = v[k][1];
      Lflip[k][1] = box - v[k][0];
    }
  }
  v = Lflip;
  r = -r;
}

int main() {

  ofstream outfile("Rij.txt");
  ofstream fout("Lambda.txt");

  L.resize(N); // N position vectors
  for (int i = 0; i < N; ++i) {
    L[i].resize(2, 0);
  } // each has 2 components
  // data.resize(iterations);

  // L = {{1., 1.}, {2., 2.}, {3., 3.}, {3., 1.}};
  for (int i = 0; i < L.size(); ++i) {
    L[i][0] = ranReal() * box;
    L[i][1] = ranReal() * box;
  }
  for (int i = 0; i < L.size(); ++i) {
    cout << L[i][0] << " , " << L[i][1] << "\n";
    cout << "---------------------------------"
         << "\n";
  }
  int rot = 1;
  std::uniform_int_distribution<int> unif3(0, L.size() - 1);
  std::function<int()> choice = std::bind(unif3, std::ref(mersenne));
  //  nextParticle.resize(numbOfNextParticles);    // #nextParticle position
  // vectors
  double ltilde = 0.9;

  for (int iter = 0; iter < iterations; ++iter) {

    int i = ran0_Nmin1();
    int j = (i + ran1_Nmin1()) % N;
    data.push_back(dist(L[i], L[j]));

    // cout << "i= " << i << " "
    //     << "j= " << j << " "
    //     << "L[i][0], L[i][1]= " << L[i][0] << ", " << L[i][1] << " "
    //     << "L[j][0], L[j][1]= " << L[j][0] << ", " << L[j][1] << " "
    //     << "dist( L[i],L[j] )= " << dist(L[i], L[j]) << endl;

    if (ranBit() < 1) {
      flipConf(L, rot);
    }
    zeroDistanceToGo = ltilde;

    int nextParticleIndex = choice();

    nextParticle[0] = L[nextParticleIndex][0];
    nextParticle[1] = L[nextParticleIndex][1];

    // here comes the while loop
    while (zeroDistanceToGo > 0.0) {

      L.erase(L.begin() + nextParticleIndex);
      // L.resize( L.size() );

      for (int i = 0; i < L.size(); ++i) {
        L[i][0] = mod(L[i][0] - nextParticle[0] + box, box);
        L[i][1] = mod(L[i][1] - nextParticle[1] + box, box);
        //        L[i][0] = mod((L[i][0] - nextParticle[0]), box);
        //        L[i][1] = mod((L[i][1] - nextParticle[1]), box);
      }

      double nextPosition = std::numeric_limits<double>::infinity();
      double a = std::numeric_limits<double>::infinity();
      nextParticle[0] = a;
      nextParticle[1] = 0.0;
      double currentPosition = 0.0;

      for (int i = 0; i < L.size(); ++i) {
        double xImage[2];
        xImage[0] = L[i][0];
        xImage[1] = L[i][1];

        if (L[i][0] > box / 2) {
          xImage[0] = L[i][0] - box;
        }

        if (L[i][1] > box / 2) {
          xImage[1] = L[i][1] - box;
        }

        if (abs(xImage[1]) < 1.0) {
          double xDummy = 0.0;
          xDummy = xImage[0] - sqrt(d - (xImage[1] * xImage[1]));
          if (xDummy > 0.0 && xDummy < min(zeroDistanceToGo, nextPosition)) {
            nextPosition = xDummy;
            nextParticle[0] = L[i][0];
            nextParticle[1] = L[i][1];
            nextParticleIndex = i;
          }
        }
      }
      double distanceToNextEvent = nextPosition - currentPosition;
      //  cout << distanceToNextEvent << "\n";
      if (zeroDistanceToGo < distanceToNextEvent) {
        currentPosition += zeroDistanceToGo;
        L.push_back({currentPosition, 0.0});
        break;
      } else {
        counts++;
        lambda0 += distanceToNextEvent;
        lambdaInterim.push_back(distanceToNextEvent);
        currentPosition += distanceToNextEvent;
        zeroDistanceToGo -= distanceToNextEvent;
        L.push_back({currentPosition, 0.0});
      }
    }
  }
  //  cout << "----------------------" << endl;
  for (int i = 0; i < data.size(); ++i) {
    //  cout /*<< std::setprecision(12) << i << "\t" */ << data[i] << "\n";
    outfile << std::setprecision(12) << data[i] << "\n";
  }
  lambda0 = lambda0 / counts;
  cout << "lambda0= " << lambda0 << endl;
  for (int i = 0; i < lambdaInterim.size(); i++) {
    lambda.push_back(lambdaInterim[i] / lambda0);
    // cout << "lambda[i]= " << " "<< lambda[i] << "\n";
  }
  for (int i = 0; i < lambda.size(); ++i) {
    fout << std::setprecision(12) << lambda[i] << "\n";
    // cout << std::setprecision(12) << lambda[i] << "\n";
  }
  return 0;
}
