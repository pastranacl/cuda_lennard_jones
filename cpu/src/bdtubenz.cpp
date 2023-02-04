#include "bdtubenz.hpp" 
#include "allvec.hpp"
#include "cpputils.hpp"
#include "rands.hpp"


#include <ctime>        // Current day
#include <sys/time.h>     // Exect time

void tubenz_bdsim (struct enzymes *enz, 
                   struct tube *tb, 
                   struct environment *env,
                   double dt,
                   double SIM_TIME,
                   int SAVING_FREQ)
{
    
    long *idum, foornd = -time(NULL);                   // Variables for randm number generator    
    
    double **gr_enz;                                    // Gradient on each enzyme
    double dr_rand;                                     // Gaussian white noise
   struct timeval t0, tf;          // OpenMP compatible calculation
    
    gr_enz = dmatrix(enz->N_enz, 3);

    
    idum = &foornd;
    
    /*********************************************************************************/
    /*                               PREPROCESSING STEPS                             */
    /*********************************************************************************/
   
   
    // Parameters for Lennard-Jones calculations
    double SIGMA_SQ; 
    SIGMA_SQ = enz->SIGMA*enz->SIGMA;  
    enz->S6 = SIGMA_SQ*SIGMA_SQ*SIGMA_SQ;
    enz->S12 = (enz->S6)*(enz->S6);
    enz->SSQ_RC = SIGMA_SQ*TWO_POW_SIXTH_SQ;
    
    
    /*********************************************************************************/
    /*                          MAIN BROWNIAN DYNAMICS CODE                          */
    /*********************************************************************************/   
    cout << "--> Main Brownian dynamics... ";
    sleep(0.5);

    double t=0;                                          // Simulation time
    int n_sims=1;                                         // Number of simulation steps
    int n_files=1;                                         // Files saved
    
    gettimeofday(&t0, NULL);
    while(t<SIM_TIME) 
    {
        
        // Gradients
        grad(enz, tb, gr_enz);
      
        // Movement
        for(int i=0; i<(enz->N_enz); i++) {    
            
            // D influences drag and random force
            dr_rand = sqrt(2*enz->D0*dt);
            dr_rand = 0; // TODO: Noise suppresed
            enz->r_enz[i][0] -= enz->D0*gr_enz[i][0]*dt/(env->kBT) + dr_rand*gasdev(idum); 
            enz->r_enz[i][1] -= enz->D0*gr_enz[i][1]*dt/(env->kBT) + dr_rand*gasdev(idum);
            enz->r_enz[i][2] -= enz->D0*gr_enz[i][2]*dt/(env->kBT) + dr_rand*gasdev(idum); 
            
            enz->r_enz[i][0] = pbc(enz->r_enz[i][0], tb->wh);
            enz->r_enz[i][2] = pbc(enz->r_enz[i][2], tb->wh);
            
        }
        
        
        // Save the data
        if(n_sims % SAVING_FREQ == 0) {
            save_file(enz->r_enz, enz->N_enz, n_files);
            n_files++;
        }
        
        t += dt;
        n_sims++;   
    }
    gettimeofday(&tf, NULL);
    
    
    double exec_time = ((tf.tv_sec  - t0.tv_sec) * 1000000u + tf.tv_usec - t0.tv_usec)/1.0e6;
    cout << "  Simulation finished" <<  endl << "  (execution time: " << exec_time << " seconds)" << endl << endl;
    
    
    // Free memory
    free_dmatrix(gr_enz, enz->N_enz);
    
}



void grad(struct enzymes *enz, struct tube *tb, double **gr_enz)
{
    
    double dx, dy, dz, dsq; 
    double r2i, r6i, r12i;
    double fLJ;
    zeros_dmatrix(gr_enz, enz->N_enz, 3);

    // GRADIENTS between enzymes
    #pragma omp parallel for private(dx,dy,dz,dsq,r2i,r6i,r12i,fLJ)
    for(int i=0; i<enz->N_enz; i++) {
        for(int j=0; j<enz->N_enz; j++) {
            
            if(i==j) continue;
            
            double xej, yej, zej;
            double xei, yei, zei;
            xei = enz->r_enz[i][0];
            yei = enz->r_enz[i][1];
            zei = enz->r_enz[i][2];
            xej = enz->r_enz[j][0];
            yej = enz->r_enz[j][1];
            zej = enz->r_enz[j][2];
            
            dx = pbc(xei-xej, tb->wh);
            dy = yei-yej;
            dz = pbc(zei-zej, tb->wh);
            dsq = dx*dx + dy*dy + dz*dz;
            //if(dsq < enz->SSQ_RC) {
                
                r2i = 1.0/dsq;
                r6i = r2i*r2i*r2i;
                r12i = r6i*r6i;
                fLJ = 24.0*(enz->EPS_EE)*(enz->S6*r6i - 2.0*enz->S12*r12i)*r2i;  
                
                #pragma omp critical
                gr_enz[i][0] += fLJ*dx;
                #pragma omp critical
                gr_enz[i][1] += fLJ*dy; 
                #pragma omp critical
                gr_enz[i][2] += fLJ*dz;
                
            //}
        }
    }


    // GRADIENT WALLS
    for(int i=0; i<enz->N_enz; i++) {
        
        // GRADIENT RIGHT wall repulsion
        dy = enz->r_enz[i][1] - 0.5*tb->L;
        dsq = dy*dy;
        if(dsq < enz->SSQ_RC) {
            r2i = 1.0/dsq;
            r6i = r2i*r2i*r2i;
            r12i = r6i*r6i;
            
            fLJ = 24.0*(enz->EPS_EWALL)*(enz->S6*r6i - 2.0*enz->S12*r12i)*r2i;               
            gr_enz[i][1] += fLJ*dy; 
        }

        // GRADIENT LEFT wall repuslion
        dy = enz->r_enz[i][1] + 0.5*tb->L;
        dsq = dy*dy;
        if(dsq < enz->SSQ_RC) {
            r2i = 1.0/dsq;
            r6i = r2i*r2i*r2i;
            r12i = r6i*r6i;
            
            fLJ = 24.0*(enz->EPS_EWALL)*(enz->S6*r6i - 2.0*enz->S12*r12i)*r2i;               
            gr_enz[i][1] += fLJ*dy; 
        }
    }
    
    
}




void save_file(double **r_enz, int N_enz, int n_files)
{
    std::string FULL_PATH_ENZYMES;                 
    
    // Path and file name structured string for enzyme coordinates
    FULL_PATH_ENZYMES = strcat(ENZYMES_PATH, int2str(n_files));
    FULL_PATH_ENZYMES = strcat(FULL_PATH_ENZYMES, ENZYMES_FNAME);
    ofstream fid_enzymes_coords;
    fid_enzymes_coords.open (FULL_PATH_ENZYMES, ofstream::out);
    
    fid_enzymes_coords << N_enz << "\n\n";
    for(int i=0; i<N_enz; i++)
        fid_enzymes_coords << r_enz[i][0] << "\t" << r_enz[i][1] << "\t" << r_enz[i][2] << "\n";
    fid_enzymes_coords.close();
    
    
}

