#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream> 
#include <algorithm>
#include <unistd.h>
#include <iomanip>          // std::setprecision




#ifndef PI
#define PI 3.141592653589793238462643383279502884197169
#endif

#define TWO_POW_SIXTH_SQ      1.259921049894873


#define ENZYMES_PATH                "./enzymes/"
#define ENZYMES_FNAME               "_enzymes.xyz"



// Standard namespace
using std::cout;
using std::cin;
using std::endl;
using std::ofstream; 
using std::ifstream; 
using std::string;
using std::getline;
using std::stringstream;
using std::nothrow; // Do not send exception but a null pointer in new




struct enzymes
{
    double **r_enz;         // Coordinates of the enzymes
    int N_enz;              // Number of enzymes (to be determined from concentration)
    double ENZ_CONC;        // Enzyme concentration on the vesicle 
    
    double K_M;             // Michaelis-Menten constant for the enzyme
    double ALPHA;           // Max factor change diffsuion constant with the substrate concentration
    double D0;          // Base diffusion constant for the enzyme
    
    double SIGMA;       // Diameter of the enzyme
    double EPS_EE;          // Depth of the repulsion potential in the enzyme-enzyme interaction
    double EPS_EWALL;        // Depth of the repulsion potential in the enzyme-membrane interaction
    
    
    double SSQ_RC;
    double S6;
    double S12;
};


struct tube
{
    double L;           // Lenght of the tube
    double wh;          // Width=Hieght of the tube
    double V;           // Volume
};


struct environment
{
    double kBT;             // Temperature of the bath [pN nm]
    double S0_MAX;          // Max concentration of the substrate [nM]
};




void tubenz_bdsim (struct enzymes *enz, 
                   struct tube *tb, 
                   struct environment *env,
                   double dt,
                   double SIM_TIME,
                   int SAVING_FREQ);


// Gradient
void grad(struct enzymes *enz, struct tube *tb, double *gr_enz);

// Save
void save_file(double **r_enz, int N_enz, int n_files);

// Periodic bounday conditions
inline double pbc(double x, double w)
{
    double xp=x;
    if(x < -w*0.5) xp +=w;
    if(x >= w*0.5) xp -=w;
    return xp;
}






