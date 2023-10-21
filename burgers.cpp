/*
Burger equation: using pseduo-spectral FG method:
∂u/∂t + u * ∂u/∂x = ν * ∂²u/∂x²
Libraries: Eigen and FFTW3
Compile: g++ -std=c++17 -o binaryfilename filename.cpp -I /path/to/Eigen/ -lfftw3
*/


#define _USE_MATH_DEFINES
#include <fftw3.h>
#include <iostream>
#include <Eigen/Dense>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <vector>
#include <fstream>

using namespace std;
using namespace Eigen;

int p = 10;
int N = pow(2,p);
int Nh = N/2.0 + 1;
double L = 2.0* M_PI;
double dx = L/N;
const double nu = 0.01;
double t_start = 0.0;
double t_end =   2.0;
double t = t_start;
double dt = 0.0001;
//Real vector
VectorXd x(N);
VectorXd ux(N);
VectorXd dux(N);
VectorXd ddux(N); 
VectorXd force_un(N);
VectorXd force_uo(N);

//complex vectors
VectorXd k(Nh); 
VectorXcd uk(Nh);
VectorXcd duk(Nh);
VectorXcd dduk(Nh);


int main() {

    ofstream outputFile;
    int outputInterval = 100;
    int outputCounter = 0;
  //initial condition
    for(int i=0;i<N;++i){
      x(i) = i*dx;
      ux(i) = sin(x(i));
    }


while(t<t_end){

fftw_plan forward_plan = fftw_plan_dft_r2c_1d(N, ux.data(), reinterpret_cast<fftw_complex*>(uk.data()), FFTW_ESTIMATE);

fftw_execute(forward_plan);
// Calculate k values,duk and dduk in spectral space
    for (int i = 0; i < Nh; ++i) {
        //k(i) = complex<double>((2.0 * M_PI * i) / L,0.0);
        k(i) = (2.0 * M_PI * i) / L;
        duk(i) = complex<double>(0.0, 1.0) * k(i) * uk(i);
        dduk(i) = -k(i) * k(i) * uk(i);
    }
    
    fftw_plan dux_plan = fftw_plan_dft_c2r_1d(N, reinterpret_cast<fftw_complex*>(duk.data()), dux.data(), FFTW_ESTIMATE);
    fftw_plan ddux_plan = fftw_plan_dft_c2r_1d(N, reinterpret_cast<fftw_complex*>(dduk.data()), ddux.data(), FFTW_ESTIMATE);
    fftw_execute(dux_plan);
    fftw_execute(ddux_plan);


 // Scale the result of the IFFT to match the original U
     dux /= N;
     ddux /=N;
    
  for (int i=0;i<N;++i){
       force_un(i) = -ux(i)*dux(i) + nu * ddux(i);
       if (t==0){
            ux(i) = ux(i) + dt * force_un(i);
       }
       else {
        ux(i) = ux(i) + dt * ((3.0/2.0) * force_un(i) - (1.0/2.0)* force_uo(i)); //adam bashforth scheme
       }
       force_uo(i) = force_un(i);
   }

  t+=dt;
    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(dux_plan);
    fftw_destroy_plan(ddux_plan);
  outputCounter++;

          // Save data every 100 time steps
        if (outputCounter == outputInterval) {
            outputCounter = 0;
            // Open a new output file with a name based on the current time step
            outputFile.open("t" + to_string(static_cast<int>(t * 1000)) + ".txt");

            // Save x(i) and ux(i) to the output file
            for (int i = 0; i < N; ++i) {
                outputFile << x(i) << " " << ux(i) << endl;
            }

            // Close the output file
            outputFile.close();
        }
    
    
      // Clean up the plans

    }
    

return 0;
}
            // Clean up the plans

    



