#include <iostream>
#include "/Users/yingqiuzhang/Documents/DSPtest/DSP/src/DspCInterface.h"
#include "/Users/yingqiuzhang/Documents/DSPtest/DSP/testexp/Example2/Randomclass.h"

#define N 30 // 30 facilities
#define M 100 //100 customers
//#define NS 10 //number of scenarios

/** This is a sample of stochastic quadratic facility location problem
 * fixed open cost, and demand cost
 * stochastic demand
*/

int main(int argc, char * argv[]){
    int seed1, seed2, NS;

    string strseed1(argv[1]);
    string strseed2(argv[2]);
    string strsce(argv[3]);
    //string strN(argv[4]);
    //string strM(argv[5]);

    //create environment variable
    DspApiEnv* env = createEnv();
    seed1 = atoi(strseed1.c_str());
    seed2 = atoi(strseed2.c_str());
    NS = atoi(strsce.c_str());

    Random R;
    R.RandomInitialise(seed1, seed2);

    double probability[NS];             // probability      
    for (int i=0; i<NS; i++){
        probability[i] = (double)1/NS;
    }
    

    /** First-stage linear cost, \sum cx */
    double c[N];
    for (int i=0; i<N; i++){
        c[i] = R.RandomDouble(100, 200);
    }

    /** First-stage variable type */
    char xtype[N];
    double xlbd[N];
    double xubd[N];
    for (int i=0; i<N; i++){
        xtype[i]='B';
        xlbd[i] = 0;
        xubd[i] = 1;
    }

    /** Second-stage variable information*/
    double yobj[M*N];
    char ytype[N*M];
    double ylbd[N*M];
    double yubd[N*M];
    for (int i=0; i<M*N; i++){
        ytype[i]='C';
        ylbd[i]=0;
        yubd[i]=1;
        yobj[i]=0;
    }
  
    /** second-stage quadratic information
     * y_ij, i is the facility, j is the customer
     * y_11, ......y_1M
     * y_21, ......y_2M
     */
    int quadyrowidx[M*N];
    int quadycolidx[M*N];
    double **quadyvalue = new double *[NS];
    for (int i=0; i<NS; i++){
        quadyvalue[i] = new double [M*N];
    }
    //double quadyvalue[NS][M*N];
    for (int i=0; i<M*N; i++){
        quadyrowidx[i]=N+i;
        quadycolidx[i]=N+i;
        //quadyvalue[i]=R.RandomDouble(20, 30);
    }
    for (int s=0; s<NS; s++){
        for (int i=0; i<M*N; i++){
            quadyvalue[s][i]=R.RandomDouble(20, 30);
        }
    }
    int numq=M*N;
    
    /** Second-stage constraint*/
    double rlbd[M*N+M];
    double rubd[M*N+M];

    CoinBigIndex start[M*N+M+1];
    int index[2*M*N+M*N];
    double value[2*M*N+M*N];
    for (int i=0; i<M*N; i++){
        start[i]=2*i;
        rubd[i]=0;
        rlbd[i]=-INFINITY;
    }
    for (int i=0; i<M; i++){
        start[M*N+i]=2*M*N+N*i;
        rlbd[M*N+i]=1;
        rubd[M*N+i]=1;
    }
    start[M*N+M]=2*M*N+M*N;

    // for (int i=0; i<M*N+N+1;i++){
    //     printf("start[%d] = %d", i, start[i]);
    // }
   
    int cnt=0;
    for (int i=0; i<N; i++){
        for (int j=0; j<M; j++){
            index[cnt]=i;
            value[cnt]=-1;
            cnt++;
            index[cnt]=N+j+i*M;
            value[cnt]=1;
            cnt++;
            
        }
    }

    for (int i=0; i<M; i++){
        for (int j=0; j<N; j++){
            index[cnt]=N+i+N*j;
            value[cnt]=1;
            cnt++;
            //cout << index[cnt] << " "  << value[cnt] <<endl;
        }
    }

    //set number of scenarios
    setNumberOfScenarios(env, NS);

    //set problem dimension
    setDimensions(env, N, 0, N*M, M*N+M);
    
    //load first-stage problem
    loadFirstStage(env, NULL, NULL, NULL, xlbd, xubd, xtype, c, NULL, NULL);
  

    //load second-stage problem corresponding to each sceanrios
    for (int i=0; i<NS; i++){
        loadQuadraticSecondStage(env, i, probability[i], start, index, value, ylbd, yubd, ytype, yobj, quadyrowidx, quadycolidx, quadyvalue[i], numq, rlbd, rubd);
    }
    for (int i=0; i<NS; i++){
        delete [] quadyvalue[i];
    }
    delete [] quadyvalue;

    setIntParam(env, "DD/MASTER/SOLVER", 4);
    setIntParam(env, "DD/SUB/SOLVER", 4);
    setIntParam(env, "DE/SOLVER", 4);
    setIntParam(env, "DD/FEAS_CUTS", 0);
    setIntParam(env, "DD/OPT_CUTS", 0);
    setIntParam(env, "DD/MASTER_ALGO", 4);
    setIntParam(env, "DD/ITER_LIM", 100);

    printf("writeMPs\n");
    //write out MPS files
    writeMps(env, "farmersQP");

    //print model
    //printModel(env);

    //set solver
    

    //solve the problem using dual decomposition
    solveDd(env);
    
    //print primal solution
    double primobj;
    primobj=getPrimalBound(env);
    //printf("optimal value = %f\n", prim);
    double solution[N+N*M*NS];
    getPrimalSolution(env, N+N*M*NS, solution);
    
    double dualobj;
    dualobj = getDualBound(env);


    // get cpu time
    double cputime;
    cputime = getCpuTime(env);

    // get wall time
    double walltime;
    walltime = getWallTime(env);

    //for (int i=0; i<N+N*M*NS;i++){
    //    printf("solution = %f\n", solution[i]);
    //}
    ofstream out;
    out.open("qUFLdd.csv", ios::app);
    //seed1 seed2 NS
    out << argv[1] << "," << argv[2] << "," << NS << "," << primobj << "," << dualobj << "," << cputime << "," << walltime << endl;
    out.close();

    //free memory
    freeSolver(env);
    freeModel(env);
    freeEnv(env);

}
