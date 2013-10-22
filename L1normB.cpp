#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "string.h"

void L1norm(double* prob1, double *prob2, mwSize dist_size, double* L1_value);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    double* prob1, *prob2;
    mwSize dist_size;
    
    prob1 = mxGetPr(prhs[0]);
    prob2 = mxGetPr(prhs[1]);
    dist_size = mxGetM(prhs[0]);
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double* L1_value = mxGetPr(plhs[0]);
    
    
    //mxArray *binaryValueArray = mxCreateDoubleMatrix(binaryValueSize, 1, mxREAL);

    L1(prob1, prob2, dist_size, L1_value);
}

void L1(double* prob1, double *prob2, mwSize dist_size, double* L1_value)
{
    
    double Ltemp = 0;
    for (int i = 0; i < dist_size; i++) {     
       
        Ltemp += fabs(prob1[i] - prob2[i]);
    }
    if (Ltemp < 1e-15){
    Ltemp = 0;
    }
    *L1_value = Ltemp;
}


// function D = L1norm(prob,prob2)
// 
// %% compute distance based on L1 norm divided by two to approximate the earth movers
// % prob,prob2 : Two distributions
// 
// if min(min(prob),min(prob2)) < 0
//     warning('Either or both distributions contain negative values');
// end
// 
// if (length(prob)~= length(prob2))
//     warning('Distributions do not have the same support');
// end
// 
// D = sum(abs(prob - prob2))/2;
// if D < 1e-15
//     D = 0;
// end