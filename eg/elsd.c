/*
 * =============================================================
 * elsd.c
 *
 * This function is for calculating
 *    E[logZ]  Z=\sum_i {w_k N_i^2} / \sum_i {N_i^2}
 *    where N_i are independent normal and w_i's are weights
 *
 *   Written by Reshad Hossein
 * =============================================================
 */

// USe freeInternalPtrs for further use
#include "mex.h"
#include <math.h>

/*
 * Input parameters prhs
 * prhs[0] : eigenvalues of the matrix
 * prhs[1] : Length of computation
 */

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    double *di, *dv, *dvn, *c, *d, b, value, v;
    int L, k, m, n;
    
    if(nrhs != 2 || nlhs > 1)
        mexErrMsgTxt("Usage: O = elsd(V,L)");
    
    /* Vector containing weights */
    di = mxGetPr(prhs[0]);
    
    /* prhs[0] is first argument.
     * mxGetPr returns double*  (data, col-major)
     * mxGetM returns int  (rows)
     * mxGetN returns int  (cols)
     */
    /* m = rows(T) */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    
    if(m != 1 && n !=1) mexErrMsgTxt("weights are vectors");
    if(m < n) m=n;
    
    /* The length for cutting the computation */
    L = mxGetScalar(prhs[1]);
    if(L <= 0) 
        mexErrMsgTxt("Second argument must be a positive integer");
    for (k = 0 ; k <= m - 1; k++) 
        if( *(di+k) <= 0 ) mexErrMsgTxt("Weights must be positive");

    /* Allocation memory for some variables */
    dv = mxCalloc(m, sizeof(double));
    dvn = mxCalloc(m, sizeof(double));
    c = mxCalloc(L, sizeof(double));
    d = mxCalloc(L, sizeof(double));
    
    // Sorting input vector
    for (k = 0 ; k <= m - 1; k++) *(dv+k) = *(di+k);
    for (k = 1 ; k <= m - 1; k++) {
        n = k;
        //mexPrintf("value %i.\n",n);
        while ( n > 0 && *(dv+n) < *(dv+n-1)) {
            value      = *(dv+n);
            *(dv+n)    = *(dv+n-1);
            *(dv+n-1)  = value;
            n--;
        }
    }
    
    // Calculating b parameter 
    b = 2 * *(dv) * *(dv+m-1) / (*(dv) + *(dv+m-1));
    
    // Calculating b/dv
    for (k = 0 ; k <= m - 1; k++) *(dv+k) = b / *(dv+k);
    
    // calculating c0
    value = sqrt( *dv );
    for (k = 1 ; k <= m - 1; k++) value *=  sqrt( *(dv+k) );
    *c = value;
    
    // Calculating d0
    for (k = 0 ; k <= m - 1; k++) *(dv+k) = 1 - *(dv+k);
    for (k = 0 ; k <= m - 1; k++) *(dvn+k) = *(dv+k);
    value = *dv;
    for (k = 1 ; k <= m - 1; k++) value +=  *(dv+k) ;
    *d = value;
    
    // main loop for calculating c and d as expressed in paper
    for (k = 1 ; k <= L - 1; k++){
        value = *c * *(d+k-1);
        for (n = 1; n <= k-1; n++)   value = value + *(c+n) * *(d+k-1-n);
        v = k; // conversion int to double otherwise 1/(2*k) = 0
        *(c+k) =  1/(2*v) * value;
        for (n = 0 ; n <= m - 1; n++) *(dvn+n) = *(dv+n) * *(dvn+n);
        value = *dvn;
        for (n = 1 ; n <= m - 1; n++) value +=  *(dvn+n) ;
        *(d+k) =  value;
    }

    // calculating cumulative c
    for (k = 1 ; k <= L - 1; k++) *(c+k) = *(c+k-1) + *(c+k);
    
    // 1-c
    value = m / 2;
    for (k = 0 ; k <= L - 1; k++) *(c+k) =  1/(value+k) * ( 1 - *(c+k) );
    
    // suming up
    value = *c;
    for (k = 1 ; k <= L - 1; k++) value = value + *(c+k);
    
    value = value + log(b);

    /* Set the output pointer to the output matrix. */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(plhs[0]) = value;
    
    /* Free alocated memory and return*/
    mxFree(dv);
    mxFree(d);
    mxFree(c);
    mxFree(dvn);
    return;
}