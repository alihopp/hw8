#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

void Hamilton(double *p1, double *p2, double *q1, double *q2, double *H, const double dt, const int N, const double e){
  ofstream out ("data.txt");

  //initial condition
  q1[0] = 1-e;
  q2[0] = 0;
  p1[0] = 0;
  p2[0] = sqrt((1+e)/(1-e));

  //0th Hamilton
  H[0] = 0.5*(p1[0]*p1[0] + p2[0]*p2[0]) - 1.0/sqrt(q1[0]*q1[0] + q2[0]*q2[0]);

  // saving 0th data
  out <<0*dt<< "\t" <<p1[0]<< "\t" <<p2[0]<< "\t" <<q1[0]<< "\t" <<q2[0]<< "\t" <<H[0]<< endl;


  //Euler-method with Hamilton and saving the data
  for (int i=0; i<N-1; i++){
    //Euler-Method
    p1[i+1] = p1[i] - dt*pow((q1[i]*q1[i])+(q2[i]*q2[i]),-3/2);
    p2[i+1] = p2[i] - dt*pow((q1[i]*q1[i])+(q2[i]*q2[i]),-3/2);

    q1[i+1] = q1[i] - dt*p1[i+1];
    q2[i+1] = q2[i] - dt*p2[i+1];
    
    //Hamilton
    H[i+1] = 0.5*(p1[i+1]*p1[i+1] + p2[i+1]*p2[i+1]) - 1.0/sqrt(q1[i+1]*q1[i+1] + q2[i+1]*q2[i+1]);
    
    // saving (N-1)th data every loop
    out <<(1+i)*dt<< "\t" <<p1[i+1]<< "\t" <<p2[i+1]<< "\t" <<q1[i+1]<< "\t" <<q2[i+1]<< "\t" <<H[i+1]<< endl;
  }

  out.close();
}

int main(void){
  const double e =0.6;  
  const double tbeg = 0.0;
  const double tend = 20*M_PI;
  const double dt = 0.0005;
  const int N =((tend-tbeg)/dt)+1;
  double p1[N], p2[N], q1[N], q2[N], H[N];

  Hamilton(p1,p2,q1,q2,H,dt,N,e);
  // data(q1,q2,p1,p2,H,dt,N);

  return 0;
}

