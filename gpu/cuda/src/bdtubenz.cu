#include "bdtubenz.hpp"
#include "allvec.hpp"
#include "cpputils.hpp"
#include "rands.hpp"


void tubenz_bdsim (struct enzymes *enz, 
                   struct tube *tb, 
                   struct environment *env,
                   double dt,
                   double SIM_TIME,
                   int SAVING_FREQ)
{
    double *gr_enz;                                    // Gradient on each enzyme
    double dr_rand;                                     // Gaussian white noise
    long *idum, foornd = -time(NULL);                   // Variables for randm number generator    
    struct timeval t0, tf;                              // OpenMP compatible calculation
    const int NUM_BLOCKS = (enz->N_enz + BLOCK_SIZE - 1) // Number of blocks for cuda
                           / BLOCK_SIZE;

    //gr_enz = dmatrix(enz->N_enz, 3);    
    gr_enz = dvector(enz->N_enz*3);    
    
    //size_t size = 3*(enz->N_enz)*sizeof(double);
    //cudaError_t error = cudaMallocManaged(&gr_enz, size);

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

    enzymes *d_enz;
    tube  *d_tb;
    double *d_gr_enz;
    double *d_r_enz;
    double *h_r_enz = enz->r_enz; // Use of intermedium variable

    cudaMalloc(&d_enz, sizeof(enzymes));
    cudaMalloc(&d_tb, sizeof(tube));
    cudaMalloc(&d_gr_enz, 3*enz->N_enz*sizeof(double));
    cudaMalloc(&d_r_enz, 3*enz->N_enz*sizeof(double));

    // We only use the structures to get the parameters
    cudaMemcpy(d_enz, enz, sizeof(enzymes), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tb, tb, sizeof(tube), cudaMemcpyHostToDevice);
    
    gettimeofday(&t0, NULL);
    while(t<SIM_TIME) 
    {
        
        // Gradients
        /* Important remark: Copy a host structure to the device including the 
            arrays contained on the structure is not trivial. I have not been able
            to read the data of the array, though the structure is copied given that
            I have acess to non-array variables. Thus, to overcome that limitation
            I look for alternatives by taking out the arrays
         See  https://stackoverflow.com/questions/30584392/how-to-pass-struct-containing-array-to-the-kernel-in-cuda
        */
        cudaMemcpy(d_r_enz, h_r_enz, 3*enz->N_enz*sizeof(double), cudaMemcpyHostToDevice);
        
        // Method 1
        //grad<<<NUM_BLOCKS, BLOCK_SIZE>>>(d_r_enz, d_enz, d_tb, d_gr_enz);

        /*  Method 2
        grad_direct<<<NUM_BLOCKS, BLOCK_SIZE>>>(d_r_enz, enz->N_enz,  
                                                enz->SSQ_RC, 
                                                enz->EPS_EE, 
                                                enz->EPS_EWALL,
                                                 enz->S6, 
                                                 enz->S12, 
                                                tb->L, tb->wh,
                                                d_gr_enz);
        */

        // Method 3. 2D array exploration HIGLY LIMITED
        dim3 dimBlock(32, 32);// A block of 32 x 32 threads (max in 1024 threads for a block)
        dim3 dimGrid(1, 1);  // Single grid
        grad_direct_2d<<<dimGrid, dimBlock>>>(d_r_enz, enz->N_enz,  
                                                enz->SSQ_RC, 
                                                enz->EPS_EE, 
                                                enz->EPS_EWALL, 
                                                enz->S6, enz->S12, 
                                                tb->L, tb->wh,
                                                d_gr_enz);
        
        // Syncrhonize and copy back
        cudaDeviceSynchronize();
        cudaMemcpy(gr_enz, d_gr_enz, 3*enz->N_enz*sizeof(double), cudaMemcpyDeviceToHost);
        


        // Movement of enzymes
        for(int i=0; i<(enz->N_enz); i++) {  
            dr_rand = sqrt(2*enz->D0*dt);
            dr_rand = 0; // TODO: NO NOISE!
            enz->r_enz[i*3+0] -= enz->D0*gr_enz[i*3 + 0]*dt/(env->kBT) + dr_rand*gasdev(idum); 
            enz->r_enz[i*3+1] -= enz->D0*gr_enz[i*3 + 1]*dt/(env->kBT) + dr_rand*gasdev(idum);
            enz->r_enz[i*3+2] -= enz->D0*gr_enz[i*3 + 2]*dt/(env->kBT) + dr_rand*gasdev(idum); 
            
            enz->r_enz[i+0] = pbc(enz->r_enz[i+0], tb->wh);
            enz->r_enz[i+2] = pbc(enz->r_enz[i+2], tb->wh);
        }
        
        
        // Save the data
        if(n_sims % SAVING_FREQ == 0) {
            save_file(enz->r_enz, enz->N_enz, n_files);
            n_files++;
        }
        
        t += dt;
        n_sims++;   
	    //break;
    }
    
    
    gettimeofday(&tf, NULL);
    double exec_time = ((tf.tv_sec  - t0.tv_sec) * 1000000u + tf.tv_usec - t0.tv_usec)/1.0e6;
    cout << "  Simulation finished" <<  endl << "  (execution time: " << exec_time << " seconds)" << endl << endl;
    
    
    // Free memory
    //free_dmatrix(gr_enz, enz->N_enz);
    // Free the managed memory
    // It'd great if you free everything you have allocated, isn't it?
    cudaFree(gr_enz);

}


void save_file(double *r_enz, int N_enz, int n_files)
{
    std::string FULL_PATH_ENZYMES;                 
    
    // Path and file name structured string for enzyme coordinates
    FULL_PATH_ENZYMES = strcat(ENZYMES_PATH, int2str(n_files));
    FULL_PATH_ENZYMES = strcat(FULL_PATH_ENZYMES, ENZYMES_FNAME);
    ofstream fid_enzymes_coords;
    fid_enzymes_coords.open (FULL_PATH_ENZYMES, ofstream::out);
    
    fid_enzymes_coords << N_enz << "\n\n";
    for(int i=0; i<N_enz; i++)
        fid_enzymes_coords << std::setprecision(9)
                           << r_enz[i*3+0] << "\t" 
                           << r_enz[i*3+1] << "\t"
                           << r_enz[i*3+2] << "\n";
                           
    fid_enzymes_coords.close();
    
}


/*************************************************************************************/
/*                                      CUDA CODES                                   */
/*************************************************************************************/



/**************************************************************************/
/*  We pass the array independently and all the structures. Note that the 
    if the structures contain arrays you have to send them to the GPU 
    memory, with the ensuing time consumption.
/**************************************************************************/
__global__ void grad(double *r_enz, enzymes *enz, tube *tb, double *gr_enz)
{

    
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= enz->N_enz) return;
    
    double dx, dy, dz, dsq; 
    double r2i, r6i, r12i;
    double fLJ;

    gr_enz[i*3+0] = 0; 
    gr_enz[i*3+1] = 0; 
    gr_enz[i*3+2] = 0; 
   
    // GRADIENTS between enzymes
    for(int j=0; j < enz->N_enz; j++) 
    {
        if(i==j) continue;
       
        double xej, yej, zej;
        double xei, yei, zei;
        xei = r_enz[i*3+0];
        yei = r_enz[i*3+1];
        zei = r_enz[i*3+2];
        xej = r_enz[j*3+0];
        yej = r_enz[j*3+1];
        zej = r_enz[j*3+2];
        
        dx = cupbc(xei-xej, tb->wh);
        dz = cupbc(zei-zej, tb->wh);
        dy = yei-yej;
        dsq = dx*dx + dy*dy + dz*dz;

        if(dsq < enz->SSQ_RC) {
            
            r2i = 1.0/dsq;
            r6i = r2i*r2i*r2i;
            r12i = r6i*r6i;
            fLJ = 24.0*(enz->EPS_EE)*(enz->S6*r6i - 2.0*enz->S12*r12i)*r2i;  
    
            gr_enz[3*i] += fLJ*dx;
            gr_enz[3*i+1] += fLJ*dy; 
            gr_enz[3*i+2] += fLJ*dz;    
	    }
    }

    // GRADIENT RIGHT wall repulsion
    dy = r_enz[i*3+1] - 0.5*tb->L;
    dsq = dy*dy;
    if(dsq < enz->SSQ_RC) {
        r2i = 1.0/dsq;
        r6i = r2i*r2i*r2i;
        r12i = r6i*r6i;
        
        fLJ = 24.0*(enz->EPS_EWALL)*(enz->S6*r6i - 2.0*enz->S12*r12i)*r2i;               
        gr_enz[3*i+1] += fLJ*dy; 
    }

    // GRADIENT LEFT wall repuslion
    dy = r_enz[i*3+1] + 0.5*tb->L;
    dsq = dy*dy;
    if(dsq < enz->SSQ_RC) {
        r2i = 1.0/dsq;
        r6i = r2i*r2i*r2i;
        r12i = r6i*r6i;
        
        fLJ = 24.0*(enz->EPS_EWALL)*(enz->S6*r6i - 2.0*enz->S12*r12i)*r2i;               
        gr_enz[3*i+1] += fLJ*dy; 
    }
    //__syncthreads();
}


/**************************************************************************/
/*  Improved version of grad. Since we require to send r_enz, we avoid 
   sending the structure, which has r_enz (but I am not able to access 
   to it), and we send the main ingredients that are necessary for the 
   function.
/**************************************************************************/
__global__ void grad_direct(double *r_enz,
                           int N_enz,
                           double enz_SSQ_RC,
                           double enz_EPS_EE,
                           double enz_EPS_EWALL,
                           double enz_S6,
                           double enz_S12, 
                           double tb_L,
                           double tb_wh,
                           double *gr_enz)
{
    

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N_enz) return;
    
    double dx, dy, dz, dsq; 
    double r2i, r6i, r12i;
    double fLJ;

    gr_enz[i*3+0] = 0; 
    gr_enz[i*3+1] = 0; 
    gr_enz[i*3+2] = 0; 
    
    // GRADIENTS between enzymes
    for(int j=0; j < N_enz; j++) 
    {
        if(i==j) continue;
       
        double xej, yej, zej;
        double xei, yei, zei;
        xei = r_enz[i*3+0];
        yei = r_enz[i*3+1];
        zei = r_enz[i*3+2];
        xej = r_enz[j*3+0];
        yej = r_enz[j*3+1];
        zej = r_enz[j*3+2];
        
        dx = cupbc(xei-xej, tb_wh);
        dz = cupbc(zei-zej, tb_wh);
        dy = yei-yej;
        dsq = dx*dx + dy*dy + dz*dz;

        if(dsq < enz_SSQ_RC) {
            
            r2i = 1.0/dsq;
            r6i = r2i*r2i*r2i;
            r12i = r6i*r6i;
            fLJ = 24.0*(enz_EPS_EE)*(enz_S6*r6i - 2.0*enz_S12*r12i)*r2i;  
    
            gr_enz[3*i] += fLJ*dx;
            gr_enz[3*i+1] += fLJ*dy; 
            gr_enz[3*i+2] += fLJ*dz;    
	    }
    }

    // GRADIENT RIGHT wall repulsion
    dy = r_enz[i*3+1] - 0.5*tb_L;
    dsq = dy*dy;
    if(dsq < enz_SSQ_RC) {
        r2i = 1.0/dsq;
        r6i = r2i*r2i*r2i;
        r12i = r6i*r6i;
        
        fLJ = 24.0*(enz_EPS_EWALL)*(enz_S6*r6i - 2.0*enz_S12*r12i)*r2i;               
        gr_enz[3*i+1] += fLJ*dy; 
    }

    // GRADIENT LEFT wall repuslion
    dy = r_enz[i*3+1] + 0.5*tb_L;
    dsq = dy*dy;
    if(dsq < enz_SSQ_RC) {
        r2i = 1.0/dsq;
        r6i = r2i*r2i*r2i;
        r12i = r6i*r6i;
        
        fLJ = 24.0*(enz_EPS_EWALL)*(enz_S6*r6i - 2.0*enz_S12*r12i)*r2i;               
        gr_enz[3*i+1] += fLJ*dy; 
    }
    //__syncthreads();
    
}


/*****************************************************************/
/* First trial to explore the array of i-j interactions considering 
   a 2d array. This a limited approximation to 32 particles. The
   reason is that we can atomicAdd only works for data in a single
   block, and is blind to what happens in other blocks. Thus, we 
   need to limit ourself to a single block, which has maximum of 
   1024 threds, 32x32 possible interactions.
/*****************************************************************/
__global__ void grad_direct_2d(double *r_enz,
                                int N_enz,
                                double enz_SSQ_RC,
                                double enz_EPS_EE,
                                double enz_EPS_EWALL,
                                double enz_S6,
                                double enz_S12, 
                                double tb_L,
                                double tb_wh,
                                double *gr_enz)
{
    

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i >= N_enz || j >= N_enz || i==j) return;
    //printf("%d", j);

    gr_enz[i*3+0] = 0; 
    gr_enz[i*3+1] = 0; 
    gr_enz[i*3+2] = 0; 
    gr_enz[j*3+0] = 0; 
    gr_enz[j*3+1] = 0; 
    gr_enz[j*3+2] = 0; 

    __syncthreads();
    
    double dx, dy, dz, dsq; 
    double r2i, r6i, r12i;
    double fLJ;

    
    // GRADIENTS between enzymes
   
    double xej, yej, zej;
    double xei, yei, zei;
    xei = r_enz[i*3+0];
    yei = r_enz[i*3+1];
    zei = r_enz[i*3+2];
    xej = r_enz[j*3+0];
    yej = r_enz[j*3+1];
    zej = r_enz[j*3+2];
    
    dx = cupbc(xei-xej, tb_wh);
    dz = cupbc(zei-zej, tb_wh);
    dy = yei-yej;
    dsq = dx*dx + dy*dy + dz*dz;

    if(dsq < enz_SSQ_RC) {
        
        r2i = 1.0/dsq;
        r6i = r2i*r2i*r2i;
        r12i = r6i*r6i;
        fLJ = 24.0*(enz_EPS_EE)*(enz_S6*r6i - 2.0*enz_S12*r12i)*r2i;  
        atomicAdd(&gr_enz[3*i], fLJ*dx);
        atomicAdd(&gr_enz[3*i+1], fLJ*dy);
        atomicAdd(&gr_enz[3*i+2], fLJ*dz);
    }

    __syncthreads();
    /*
    if(j==1) { // This is to avoid double calculation
        // GRADIENT RIGHT wall repulsion
        dy = r_enz[i*3+1] - 0.5*tb_L;
        dsq = dy*dy;
        if(dsq < enz_SSQ_RC) {
            r2i = 1.0/dsq;
            r6i = r2i*r2i*r2i;
            r12i = r6i*r6i;
            
            fLJ = 24.0*(enz_EPS_EWALL)*(enz_S6*r6i - 2.0*enz_S12*r12i)*r2i;               
            gr_enz[3*i+1] += fLJ*dy; 
        }

        // GRADIENT LEFT wall repuslion
        dy = r_enz[i*3+1] + 0.5*tb_L;
        dsq = dy*dy;
        if(dsq < enz_SSQ_RC) {
            r2i = 1.0/dsq;
            r6i = r2i*r2i*r2i;
            r12i = r6i*r6i;
            
            fLJ = 24.0*(enz_EPS_EWALL)*(enz_S6*r6i - 2.0*enz_S12*r12i)*r2i;               
            gr_enz[3*i+1] += fLJ*dy; 
        }
    }
    */
}

