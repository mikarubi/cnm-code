#include "mex.h"                        /* for mex                  */
#include "string.h"                     /* for memcpy               */
#include "time.h"                       /* for rand initialization  */
#include "math.h"                       /* for sqrt and fabs        */

/*declare rand_mt functions*/
void init_genrand(unsigned long s);
double genrand_real2(void);

void theloop(
        double *const W[],              double  *const B[],
        const double *const D,          const int32_T *const M,
        const int32_T objS,             const int32_T objL,             const int32_T objR,
        const int32_T *const Ns,        const int32_T *const Na,
        const int32_T *const Ms,        const int32_T *const INa,       const int32_T n_asym,
        const int32_T *const I,         const int32_T *const J,
        const int32_T *const IdxN,      const int32_T *const IdxM,      const int32_T binary_cn,
        double temp,                    const double temp_decy,         const double thr1,
        const int32_T n,                const int32_T n_edge,           const int32_T m,
        const int32_T n_obj,            const int32_T time_totl,
        const int32_T time_decy,        const int32_T rng_seed,
        double *const W0[],           	double *const B0[],
        double *const HistoryMeanDelt,  double *const HistoryThrshold)
{
    int t, h, i, j, e1, e2, a, b, c, d, n_vald_node, n_vald_modu, p = binary_cn + 1;
    double  DeltS[p][2*n], DeltM[p][m*m], DeltL[p][2*n], ScleS[p], ScleM[p], ScleL[p],
            wcd_minus_wab_nrms[p], wcd_minus_wab_nrmm[p],
            rcd_minus_rab_nrml[p], rab_minus_rcd_nrml[p],
            wab, wcd, bab, bcd, summ, delt0, delt1, delti;
    
    int a1, b1, c1, d1;
    double  wc1d1_minus_wa1b1_nrms[p], wc1d1_minus_wa1b1_nrmm[p],
            rc1d1_minus_ra1b1_nrml[p], ra1b1_minus_rc1d1_nrml[p],
            wa1b1, wc1d1, ba1b1, bc1d1;
    
    init_genrand((rng_seed*100)+(time(NULL)/100));              /* seed combines time and iteration */
    
    /* initialize arrays */
    delt0 = 0;
    for(h=0; h<p; h++){
        ScleS[h] = 0;
        ScleM[h] = 0;
        ScleL[h] = 0;
        for(i=0; i<2*n; i++) DeltS[h][i] = 0;
        for(i=0; i<m*m; i++) DeltM[h][i] = 0;
        for(i=0; i<2*n; i++) DeltL[h][i] = 0;
    }
    
    /* binarize one of the matrices */
    for(i=0; i<n; i++)
        for(j=0; j<n; j++){
            W0[1][i + j*n] = (W0[1][i + j*n]>0);
            B0[1][i + j*n] = (B0[1][i + j*n]>0);
            
            W[1][i + j*n] = (W[1][i + j*n]>0);
            B[1][i + j*n] = (B[1][i + j*n]>0);
        }
    
    /* get n_vald_node */
    n_vald_node = 0;
    for(i=0; i<2*n; i++)
        if(IdxN[i])
            n_vald_node += 1;
    
    /* get n_vald_modu */
    n_vald_modu=0;
    for(i=0; i<m*m; i++)
        if(IdxM[i])
            n_vald_modu += 1;
    
    for(h=0; h<p; h++){
        if(objS){
            /* get ScleS, which depends on n_vald_node */
            for(i=0; i<n; i++){
                if(IdxN[  i])
                    for(j=0; j<n; j++)
                        ScleS[h] += W[h][i + j*n]/n_vald_node;
                if(IdxN[n+i])
                    for(j=0; j<n; j++)
                        ScleS[h] += W[h][j + i*n]/n_vald_node;
            }
            
            /* get DeltS, which depends on ScleS */
            for(i=0; i<n; i++){
                if(IdxN[  i])
                    for(j=0; j<n; j++)
                        DeltS[h][  i] += (W0[h][i + j*n] - W[h][i + j*n])/ScleS[h];
                if(IdxN[n+i])
                    for(j=0; j<n; j++)
                        DeltS[h][n+i] += (W0[h][j + i*n] - W[h][j + i*n])/ScleS[h];
            }
            
            /* cumulate delt0 */
            for(delti=0, i=0; i<2*n; i++)
                if(IdxN[i])
                    delti += fabs(DeltS[h][i])/n_vald_node;
            delt0 += delti*delti;
        }
        
        if(objL){
            /* get ScleM, which depends on n_vald_modu */
            for(i=0; i<n; i++)
                for(j=0; j<n; j++)
                    if((M[i]>=0) && (M[j]>=0) && IdxM[M[i] + M[j]*m])
                        ScleM[h] += W[h][i + j*n]/n_vald_modu;
            
            /* get DeltM, which depends on ScleM */
            for(i=0; i<n; i++)
                for(j=0; j<n; j++)
                    if((M[i]>=0) && (M[j]>=0) && IdxM[M[i] + M[j]*m])
                        DeltM[h][M[i] + M[j]*m] += (W0[h][i + j*n] - W[h][i + j*n])/ScleM[h];
            
            /* cumulate delt0 */
            for(delti=0, i=0; i<m*m; i++)
                if(IdxM[i])
                    delti += fabs(DeltM[h][i])/n_vald_modu;
            delt0 += delti*delti;
        }
        
        if(objR){
            /* get ScleL, which depends on n_vald_node */
            for(i=0; i<n; i++){
                if(IdxN[  i])
                    for(j=0; j<n; j++)
                        ScleL[h] += B[h][i + j*n]*D[i + j*n]/n_vald_node;
                if(IdxN[n+i])
                    for(j=0; j<n; j++)
                        ScleL[h] += B[h][j + i*n]*D[j + i*n]/n_vald_node;
            }
            
            for(i=0; i<n; i++){
                if(IdxN[  i])
                    for(j=0; j<n; j++)
                        DeltL[h][  i] += (B0[h][i + j*n] - B[h][i + j*n])*D[i + j*n]/ScleL[h];
                if(IdxN[n+i])
                    for(j=0; j<n; j++)
                        DeltL[h][n+i] += (B0[h][j + i*n] - B[h][j + i*n])*D[j + i*n]/ScleL[h];
            }
            
            /* cumulate delt0 */
            for(delti=0, i=0; i<2*n; i++)
                if(IdxN[i])
                    delti += fabs(DeltL[h][i])/n_vald_node;
            delt0 += delti*delti;
        }
    }
    
    delt0 = sqrt(delt0/(p*n_obj));
    
    mexPrintf("Initial delta:\t%f\n", delt0);
    if(time_totl==0)
        return;
    
    for(t=0; t<time_totl; t++){
        while(true){
            e1 = (int) (n_edge*genrand_real2());
            e2 = (int) (n_edge*genrand_real2());
            
            a  = I[e1];
            b  = J[e1];
            c  = I[e2];
            d  = J[e2];
            
            if(fabs(W0[0][a + b*n] - W0[0][c + d*n]) > 1.0e-10){
                if((Na[a] && Na[b]) == (Na[c] && Na[d]))        /* if a&b are both assymmetric     */
                    break;                                      /* c&d must be too, and vice versa */
            }
        }
        
        a1 = Ns[a] ? Ms[a]: a;                                  /* map to symmetric nodes          */
        b1 = Ns[b] ? Ms[b]: b;
        c1 = Ns[c] ? Ms[c]: c;
        d1 = Ns[d] ? Ms[d]: d;
        
        while(((a1==a) && (b1==b)) || ((c1==c) && (d1==d)) || (a1==b1) || (c1==d1)){
            if(Na[a] && Na[b]){
                a1 = INa[(int) (n_asym*genrand_real2())];
                b1 = INa[(int) (n_asym*genrand_real2())];
            }
            if(Na[c] && Na[d]){
                c1 = INa[(int) (n_asym*genrand_real2())];
                d1 = INa[(int) (n_asym*genrand_real2())];
            }
        }
        
        if( ((a1==c)&&(b1==d)) || ((c1==a)&&(d1==b)))
            continue;
        
        delt1 = 0;
        for(h=0; h<p; h++){
            if(objS){
                wcd_minus_wab_nrms[h]     = (W0[h][c + d*n] - W0[h][a + b*n])/ScleS[h];
                wc1d1_minus_wa1b1_nrms[h] = (W0[h][c1+d1*n] - W0[h][a1+b1*n])/ScleS[h];
                
                for(delti=0, i=0; i<2*n; i++){
                    if(IdxN[i]){
                        summ = 0;
                        if(i==a  || i==n+b )
                            summ += wcd_minus_wab_nrms[h];
                        if(i==c  || i==n+d )
                            summ -= wcd_minus_wab_nrms[h];
                        if(i==a1 || i==n+b1)
                            summ += wc1d1_minus_wa1b1_nrms[h];
                        if(i==c1 || i==n+d1)
                            summ -= wc1d1_minus_wa1b1_nrms[h];
                        
                        delti += fabs(DeltS[h][i]+summ)/n_vald_node;
                    }
                }
                delt1 += delti*delti;
            }
            
            if(objL){
                wcd_minus_wab_nrmm[h]     = (W0[h][c + d*n] - W0[h][a + b*n])/ScleM[h];
                wc1d1_minus_wa1b1_nrmm[h] = (W0[h][c1+d1*n] - W0[h][a1+b1*n])/ScleM[h];
                
                for(delti=0, i=0; i<m*m; i++){
                    if(IdxM[i]){
                        summ = 0;
                        if((i==M[ a]+M[ b]*m) && (M[ a]>=0) && (M[ b]>=0))
                            summ += wcd_minus_wab_nrmm[h];
                        if((i==M[ c]+M[ d]*m) && (M[ c]>=0) && (M[ d]>=0))
                            summ -= wcd_minus_wab_nrmm[h];
                        if((i==M[a1]+M[b1]*m) && (M[a1]>=0) && (M[b1]>=0))
                            summ += wc1d1_minus_wa1b1_nrmm[h];
                        if((i==M[c1]+M[d1]*m) && (M[c1]>=0) && (M[d1]>=0))
                            summ -= wc1d1_minus_wa1b1_nrmm[h];
                        
                        delti += fabs(DeltM[h][i]+summ)/n_vald_modu;
                    }
                }
                delt1 += delti*delti;
            }
            
            if(objR){
                rcd_minus_rab_nrml[h]     = (B0[h][c + d*n] - B0[h][a + b*n])*D[a + b*n]/ScleL[h];
                rab_minus_rcd_nrml[h]     = (B0[h][a + b*n] - B0[h][c + d*n])*D[c + d*n]/ScleL[h];
                rc1d1_minus_ra1b1_nrml[h] = (B0[h][c1+d1*n] - B0[h][a1+b1*n])*D[a1+b1*n]/ScleL[h];
                ra1b1_minus_rc1d1_nrml[h] = (B0[h][a1+b1*n] - B0[h][c1+d1*n])*D[c1+d1*n]/ScleL[h];
                
                for(delti=0, i=0; i<2*n; i++){
                    if(IdxN[i]){
                        summ = 0;
                        if(i==a  || i==n+b )
                            summ += rcd_minus_rab_nrml[h];
                        if(i==c  || i==n+d )
                            summ += rab_minus_rcd_nrml[h];
                        if(i==a1 || i==n+b1)
                            summ += rc1d1_minus_ra1b1_nrml[h];
                        if(i==c1 || i==n+d1)
                            summ += ra1b1_minus_rc1d1_nrml[h];
                        
                        delti += fabs(DeltL[h][i]+summ)/n_vald_node;
                    }
                }
                delt1 += delti*delti;
            }
        }
        delt1 = sqrt(delt1/(p*n_obj));
        
        if((exp((delt0 - delt1)/temp) > genrand_real2()) || (delt1 < thr1)){
            delt0 = delt1;
            for(h=0; h<p; h++){
                wab   = W0[h][a + b*n];
                wcd   = W0[h][c + d*n];
                bab   = B0[h][a + b*n];
                bcd   = B0[h][c + d*n];
                wa1b1 = W0[h][a1+b1*n];
                wc1d1 = W0[h][c1+d1*n];
                ba1b1 = B0[h][a1+b1*n];
                bc1d1 = B0[h][c1+d1*n];
                
                W0[h][a + b*n] = wcd;
                W0[h][c + d*n] = wab;
                B0[h][a + b*n] = bcd;
                B0[h][c + d*n] = bab;
                W0[h][a1+b1*n] = wc1d1;
                W0[h][c1+d1*n] = wa1b1;
                B0[h][a1+b1*n] = bc1d1;
                B0[h][c1+d1*n] = ba1b1;
                
                if(objS){
                    DeltS[h][   a] += wcd_minus_wab_nrms[h];
                    DeltS[h][n+ b] += wcd_minus_wab_nrms[h];
                    DeltS[h][   c] -= wcd_minus_wab_nrms[h];
                    DeltS[h][n+ d] -= wcd_minus_wab_nrms[h];
                    DeltS[h][  a1] += wc1d1_minus_wa1b1_nrms[h];
                    DeltS[h][n+b1] += wc1d1_minus_wa1b1_nrms[h];
                    DeltS[h][  c1] -= wc1d1_minus_wa1b1_nrms[h];
                    DeltS[h][n+d1] -= wc1d1_minus_wa1b1_nrms[h];
                }
                if(objL){
                    if((M[ a]>=0) && (M[ b]>=0))
                        DeltM[h][M[ a]+M[ b]*m] += wcd_minus_wab_nrmm[h];
                    if((M[ c]>=0) && (M[ d]>=0))
                        DeltM[h][M[ c]+M[ d]*m] -= wcd_minus_wab_nrmm[h];
                    if((M[a1]>=0) && (M[b1]>=0))
                        DeltM[h][M[a1]+M[b1]*m] += wc1d1_minus_wa1b1_nrmm[h];
                    if((M[c1]>=0) && (M[d1]>=0))
                        DeltM[h][M[c1]+M[d1]*m] -= wc1d1_minus_wa1b1_nrmm[h];
                }
                if(objR){
                    DeltL[h][   a] += rcd_minus_rab_nrml[h];
                    DeltL[h][n+ b] += rcd_minus_rab_nrml[h];
                    DeltL[h][   c] += rab_minus_rcd_nrml[h];
                    DeltL[h][n+ d] += rab_minus_rcd_nrml[h];
                    DeltL[h][  a1] += rc1d1_minus_ra1b1_nrml[h];
                    DeltL[h][n+b1] += rc1d1_minus_ra1b1_nrml[h];
                    DeltL[h][  c1] += ra1b1_minus_rc1d1_nrml[h];
                    DeltL[h][n+d1] += ra1b1_minus_rc1d1_nrml[h];
                }
            }
        }
        
        if( (t % time_decy) == 0){
            if(t > 0)
                temp *= temp_decy;
            
            HistoryMeanDelt[t/time_decy] = delt0;
            HistoryThrshold[t/time_decy] = temp;
        }
    }
    
    mexPrintf("Final delta:\t%f\n", delt0);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxArray
            *mxW0                	=  mexGetVariable("caller","W0"),
            *mxB0                 	=  mexGetVariable("caller","B0"),
            *mxW                    =  mexGetVariable("caller","W"),
            *mxB                    =  mexGetVariable("caller","B"),
            *mxW0bin             	=  mexGetVariable("caller","W0"),
            *mxB0bin              	=  mexGetVariable("caller","B0"),
            *mxWbin                 =  mexGetVariable("caller","W"),
            *mxBbin                 =  mexGetVariable("caller","B"),
            *mxHistoryMeanDelt      =  mexGetVariable("caller","HistoryMeanDelt"),
            *mxHistoryThrshold      =  mexGetVariable("caller","HistoryThreshold");
    
    double
            *const W0[]             =  {mxGetPr(mxW0),  mxGetPr(mxW0bin)},
            *const B0[]             =  {mxGetPr(mxB0),  mxGetPr(mxB0bin)},
            *const W[]              =  {mxGetPr(mxW),   mxGetPr(mxWbin)},
            *const B[]              =  {mxGetPr(mxB),   mxGetPr(mxBbin)},
            *const HistoryMeanDelt  =  mxGetPr(mxHistoryMeanDelt),
            *const HistoryThrshold  =  mxGetPr(mxHistoryThrshold);
                    
    const double
            *const D                =  mxGetPr(mexGetVariablePtr("caller","D")),
            temp                    = *mxGetPr(mexGetVariablePtr("caller","temp")),
            temp_decy               = *mxGetPr(mexGetVariablePtr("caller","temp_decy")),
            thr1                    = *mxGetPr(mexGetVariablePtr("caller","thr1"));

    const int32_T
            *const IdxN             =  (int32_T *)mxGetData(mexGetVariablePtr("caller","IdxN")),
            *const IdxM             =  (int32_T *)mxGetData(mexGetVariablePtr("caller","IdxM")),
            *const M                =  (int32_T *)mxGetData(mexGetVariablePtr("caller","M_c")),
            *const Ns               =  (int32_T *)mxGetData(mexGetVariablePtr("caller","Ns")),
            *const Na               =  (int32_T *)mxGetData(mexGetVariablePtr("caller","Na")),
            *const Ms               =  (int32_T *)mxGetData(mexGetVariablePtr("caller","Ms_c")),
            *const INa              =  (int32_T *)mxGetData(mexGetVariablePtr("caller","INa_c")),
            n_asym                  = *(int32_T *)mxGetData(mexGetVariablePtr("caller","n_asym")),
            *const I                =  (int32_T *)mxGetData(mexGetVariablePtr("caller","I_c")),
            *const J                =  (int32_T *)mxGetData(mexGetVariablePtr("caller","J_c")),
            binary_cn               = *(int32_T *)mxGetData(mexGetVariablePtr("caller","binary_cn")),
            objS                    = *(int32_T *)mxGetData(mexGetVariablePtr("caller","objS")),
            objL                    = *(int32_T *)mxGetData(mexGetVariablePtr("caller","objL")),
            objR                    = *(int32_T *)mxGetData(mexGetVariablePtr("caller","objR")),
            n                       = *(int32_T *)mxGetData(mexGetVariablePtr("caller","n")),
            m                       = *(int32_T *)mxGetData(mexGetVariablePtr("caller","m")),
            n_edge                  = *(int32_T *)mxGetData(mexGetVariablePtr("caller","n_edge")),
            n_obj                   = *(int32_T *)mxGetData(mexGetVariablePtr("caller","n_obj")),
            time_totl             	= *(int32_T *)mxGetData(mexGetVariablePtr("caller","time_totl")),
            time_decy             	= *(int32_T *)mxGetData(mexGetVariablePtr("caller","time_decy")),
            rng_seed                = *(int32_T *)mxGetData(mexGetVariablePtr("caller","rng_seed"));

    theloop(W, B, D, M, objS, objL, objR, Ns, Na, Ms, INa, n_asym,
            I, J, IdxN, IdxM, binary_cn, temp, temp_decy, thr1, n, n_edge, m, n_obj,
            time_totl, time_decy, rng_seed, W0, B0, HistoryMeanDelt, HistoryThrshold);

    mexPutVariable("caller", "W0",          mxW0);
    mexPutVariable("caller", "B0",          mxB0);
    mexPutVariable("caller", "HistoryMeanDelt",    mxHistoryMeanDelt);
    mexPutVariable("caller", "HistoryThreshold",   mxHistoryThrshold);
}
