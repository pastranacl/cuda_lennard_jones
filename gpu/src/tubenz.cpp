/****************************************************************************

    tubenz
    Copyright (C) 2021  Cesar L. Pastrana

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

****************************************************************************/

#include "tubenz.hpp"
#include "bdtubenz.hpp"
#include "allvec.hpp"
#include "rands.hpp"
#include "cpputils.hpp"


int main()
{
    
    double dt;                              // Simulation time step
    double SIM_TIME;                        // Max simulation time
    double SAVING_FREQ;                     // Saving frequency (inverse files per simulation step) 
    
    struct enzymes enz;                     // Struct of the enzyme properties
    struct tube tb;                         // Struct of the tube
    struct environment env;                 // Struct for the properties of the experiment

    
    //************************************************************************************//
    // 1. Load data                                                                       //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

    // 1.1  Load the file with the initial coordiantes of the particles
    dvector_vl r_enz;                       // Coordinates xyz for the enzyme(from file)       
    ifstream infile_coords(INIT_COORDS_FNAME);
    string line;
    
    bool errloadenz = 0;
    int l=0;                                 
    if(infile_coords.is_open()) {
        errloadenz = 0;
        while (getline(infile_coords, line))
        {
            double val;
            stringstream ss(line);
            r_enz.push_back(vector<double>(0));

            while (ss >> val)
                r_enz[l].push_back(val);
            
            ++l;
        }
        infile_coords.close();
        cout << "Enzyme coordinates found"; 
    } else {
        errloadenz = 1;
        cout << "Enzyme coordinates not found. A new configuration will be generated";
    }
    
    
    enz.N_enz = r_enz.size();
    
    
    
    // 1.3 Load parameters from file
    dt = 1.0;
    SIM_TIME = 1.0e6;  
    SAVING_FREQ = 1000;
    
    
    // std::ifstream is RAII, i.e. no need to call close
    ifstream infile_params(PARAMETERS_FILE);
    if (infile_params.is_open())
    {
        while (getline(infile_params, line)) 
        {
            line.erase(remove_if(line.begin(), line.end(), ::isspace),line.end());
            
            // Commments with the sharp symbol
            if (line[0] == '#' || line.empty()) continue; 

            auto delim_pos = line.find("=");
            auto tag = line.substr(0, delim_pos);
            auto val = line.substr(delim_pos + 1);

            // Tube properties
            if (tag == "L") tb.L = std::stod(val);      
            else if (tag == "wh") tb.wh = std::stod(val);

            
            // Enzyme properties
            else if (tag == "ENZ_CONZ") enz.ENZ_CONC = std::stod(val);
            else if (tag == "D0") enz.D0 = std::stod(val);
            else if (tag == "SIGMA") enz.SIGMA = std::stod(val);
            else if (tag == "EPS_EE") enz.EPS_EE = std::stod(val);
            else if (tag == "EPS_EWALL") enz.EPS_EWALL = std::stod(val);
            
            // Simulation properties
            else if (tag == "dt") dt = std::stod(val);
            else if (tag == "SIM_TIME") SIM_TIME = std::stod(val);
            else if (tag == "SAVING_FREQ") SAVING_FREQ = std::stod(val);
             
            // Environment properties
            else if (tag == "kBT") env.kBT = std::stod(val);
        }
        
        // In principle not necessary beucase ifstream self cleans
        infile_params.close();
    }
    else 
    {
        cout << "Error while opening the parameters file " << PARAMETERS_FILE << ". Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
    
 

    // -------------------------------------------------------------------------------
   
    
    //************************************************************************************//
    // 2. Print GUI                                                                       //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    cout << endl << endl;
    cout << "===============================================================" << endl;
    cout << "                 SIMULATION ENZYMES IN TUBE                    " << endl;
    cout << "===============================================================" << endl << endl;
    cout << "> TUBE PROPERTIES ------------------------------------------" << endl;
    cout << " - Tube length (L) = " << tb.L << " um" << endl;
    cout << " - Width/Heigth = " << tb.wh << " um" << endl;
    cout << "---------------------------------------------------------------" << endl << endl;
    
    cout << "> ENZYME PROPERTIES -------------------------------------------" << endl;
    cout << " - [enzymes] = " << enz.ENZ_CONC << " nM"  << endl;
    cout << " - Basal diffusion constant (D0) = " << enz.D0 << " um^2/s" << endl;
    cout << " - Enzyme diameter (sigma) = " << enz.SIGMA << " nm" << endl;
    cout << " - Repulsion potential depth enzyme-enzyme (EPS_EE) = " << enz.EPS_EE << " pN nm" << endl;
    cout << " - Repulsion potential depth enzyme-walls (EPS_WALL) = " << enz.EPS_EWALL << " pN nm" << endl << endl;
    
    cout << "> ENVIRONMENT PROPERTIES --------------------------------------" << endl;
    cout << " - Thermal energy (kBT) = " << env.kBT << " pN nm" << endl;
    cout << "---------------------------------------------------------------" << endl << endl;

    
    // ------------------------------------------------------------------------------------ 
    
    
    
    //************************************************************************************//
    // 3. Rescale and calculate associated quantities                                     //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    
    // The user provides the units in the more convenient form
    // The units of the simulation are
    // Length = nm;     Force = pN;  Time = us;
    
    tb.L *= 1.0e3;  // From um to nm
    tb.wh *= 1.0e3;  // From um to nm
    SIM_TIME *= 1.0e6; // From s to us

    
    // Calculate the number of enzymes, either use loaded file or position them
    tb.V =  tb.L*(tb.wh*tb.wh); // Volume in nm3
    
    
    
    if(errloadenz==1) {
        
        enz.N_enz = enz.ENZ_CONC*tb.V*AVOGADRO_PREFACTOR*1.0e-10; 
        enz.r_enz = dmatrix(enz.N_enz, 3);
        enzyme_positions(&enz, &tb);
        
        
        //FOR TESTS
        /*
        enz.N_enz = 5; 
        enz.r_enz = dmatrix(enz.N_enz, 3);
        for(int i=0;i<enz.N_enz; i++){
            enz.r_enz[i][0] = (0.5*i)*enz.SIGMA;
            enz.r_enz[i][1] = (0.25*i)*enz.SIGMA;
            enz.r_enz[i][2] = (0.0*i)*enz.SIGMA;
        }
        */

    } else {
        
        // Copy the enzymes from a dvector to "new" matrix
        enz.r_enz = dmatrix(enz.N_enz, 3);
        for(int i = 0; i<enz.N_enz; i++) {
            enz.r_enz [i][0] = r_enz[i][0];
            enz.r_enz [i][1] = r_enz[i][1];
            enz.r_enz [i][2] = r_enz[i][2];
        }
        r_enz.clear();
    }
    
    
    
    // Save initial coordinates
    ofstream fid_enzymes_coords;
    fid_enzymes_coords.open (INIT_COORDS_SAVE, ofstream::out);
    fid_enzymes_coords << enz.N_enz << "\n\n";
    for(int i=0; i<enz.N_enz; i++)
        fid_enzymes_coords << enz.r_enz[i][0] << "\t" << enz.r_enz[i][1] << "\t" << enz.r_enz[i][2] << "\n";
    fid_enzymes_coords.close();
    
    
    std::string FULL_PATH_ENZYMES;
    FULL_PATH_ENZYMES = strcat(ENZYMES_PATH, int2str(0));
    FULL_PATH_ENZYMES = strcat(FULL_PATH_ENZYMES, ENZYMES_FNAME);
    fid_enzymes_coords.open (FULL_PATH_ENZYMES, ofstream::out);
    fid_enzymes_coords << enz.N_enz << "\n\n";
    for(int i=0; i<enz.N_enz; i++)
        fid_enzymes_coords << enz.r_enz[i][0] << "\t" << enz.r_enz[i][1] << "\t" << enz.r_enz[i][2] << "\n";
    fid_enzymes_coords.close();
    
    
    
    // -------------------------------------------------------------------------------
    
    
    //************************************************************************************//
    // 4. Calling to the main simulation routine                                          //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    
    cout << "> Simulation running (#enzymes = " << enz.N_enz << "): " << endl;
    sleep(0.5);
    double exec_time;
    time_t start, end;   

    time(&start); 
    
    tubenz_bdsim(&enz, &tb, &env, dt, SIM_TIME, SAVING_FREQ);
    
    time(&end); 
    exec_time = double(end-start);
    
    
    cout << "Done" <<  endl << "  (execution time total: " << exec_time << " seconds)" << endl << endl;
    
    
    
    
    // -------------------------------------------------------------------------------------
    

    
    return 0;
}




/*****************************************************************************************/
/*                                  enzyme_positions                                     */
/* Random positioning of the enzymes inside the tube taking care of not overlaping        */
/*****************************************************************************************/

void enzyme_positions(struct enzymes *enz, struct tube *tb)
{
    long *idum, foornd = -time(NULL);        // Variables for randm number generator    
    
    const double SP=3.0;                    // 3xsigma separation
    double d_min;                           // Minimum separation between a couple of enzymes
    double maxx, minx, maxL, minL;          // Boundaries of the tube
    double sq_sigma;                        // Sigma square
        
    idum = &foornd;
    
    
    maxx  =  tb->wh/2.0 - SP*enz->SIGMA;
    minx  = -tb->wh/2.0 + SP*enz->SIGMA;
    
    maxL =  tb->L/2.0 - SP*enz->SIGMA;
    minL = -tb->L/2.0 + SP*enz->SIGMA;

    sq_sigma = (SP*enz->SIGMA)*(SP*enz->SIGMA);
    
   
    // First particle random
    enz->r_enz[0][0] = (maxx-minx)*ran1(idum) + minx;
    enz->r_enz[0][1] = (maxL-minL)*ran1(idum) + minL;
    enz->r_enz[0][2] = (maxx-minx)*ran1(idum) + minx;
    
    for(int i=1; i<enz->N_enz; i++) {
        
        enz->r_enz[i][0] = (maxx-minx)*ran1(idum) + minx;
        enz->r_enz[i][1] = (maxL-minL)*ran1(idum) + minL;
        enz->r_enz[i][2] = (maxx-minx)*ran1(idum) + minx;
        
        d_min=det_dmin(enz, i);
        
        while(d_min<sq_sigma) {
      
            enz->r_enz[i][0] = (maxx-minx)*ran1(idum) + minx;
            enz->r_enz[i][1] = (maxL-minL)*ran1(idum) + minL;
            enz->r_enz[i][2] = (maxx-minx)*ran1(idum) + minx;
    
            // Check that two enzymes are not too near
            d_min=det_dmin(enz, i);
        } 
        
    }
    
}


double det_dmin(struct enzymes *enz, int id)
{
    
    double d_min;
    double dx,dy,dz,d;
    
    dx = enz->r_enz[id][0] - enz->r_enz[0][0];
    dy = enz->r_enz[id][1] - enz->r_enz[0][1];
    dz = enz->r_enz[id][2] - enz->r_enz[0][2];
    d_min = dx*dx + dy*dy + dz*dz;
    
    
    for(int j=0; j<id; j++){
        dx = enz->r_enz[id][0] - enz->r_enz[j][0];
        dy = enz->r_enz[id][1] - enz->r_enz[j][1];
        dz = enz->r_enz[id][2] - enz->r_enz[j][2];
        d = dx*dx + dy*dy + dz*dz;
        if(d<d_min)
            d_min = d;
    }
    
    return d_min;
    
}



    

