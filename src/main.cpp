#include "../include/read_formula.h"
#include "../include/bdd_gradient.h"
#include "../include/optimize.h"
#include "../include/main.h"
#include <chrono>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <stdint.h>
#include <errno.h>
#include <limits.h>

void print_usage(char* prog_name, Parameter *param)
{
    printf( "%s [OPTIONS] InputFile (.cnf/.opb/.wcnf)\n", prog_name); 
    printf( "OPTIONS:\n");
    printf( "\t-o [INT]: Number of constriants that can be violated. (default %d)\n", param->n_violance_tolerance);
    printf( "\t-maxsat: Solve maxsat instances. (default: solve SAT problems)\n");
    printf( "\t-maxsat -wcnf: Solve weighted maxsat instances.\n");
    printf( "\t-n [INT]: Number of cores, 4 cores or more recommended. (default %d)\n", param->ncores);
    printf( "\t-s [STRING]: Continuous optimizer ('SLSQP','CG','BFGS','MMA'). (default: A portfolio of optimizers when using multicore) \n");
    printf( "\t-m [DOUBLE]: Multiplier r for adaptive constraint weight. (default %f)\n", param->multiplier);
    printf( "\t-t [INT]: Maximum number of trials T per restart for adaptive constraint weight. (default %d)\n", param->max_trial_per_start);
}

void get_parameter(int argc, char **argv, Parameter *param)
{
    Parameter _param = {
        "PORTFOLIO", // solver
         0, // k
	 1e-3, //eps
         1000, // maxiter
         2, // multiplier
         1, // verbose
        NULL, //fin
        NULL, // fin name
        0, //beta
        0, // adapt
        1, // ncores
        8, //max trial per start
        0, //L2
        0, // violance tolerance
        0, // is maxsat
        0, // is wcnf
        0 //wcnf_weight
    };

    if(argc <= 1){
        print_usage(argv[0], &_param);
        exit(0);
    }

    char **p = argv+1;
    int i;
    for(i=1; i<argc; i++, p++){
        if(!strcmp(*p, "-s")){
            if(i+1 >= argc) break;
            _param.solver = p[1];
            i++, p++;
        }else if(!strcmp(*p, "-l")){
            if(i+1  >= argc) break;
            int ret = sscanf(p[1], "%f", &_param.L2);
            if(ret != 1){
                break;
            }
            i++, p++;
        }else if(!strcmp(*p, "-e")){
            if(i+1 >= argc) break; 
            int ret = sscanf(p[1], "%lf", &_param.eps);
            if(ret != 1) break; 
            i++, p++;
        }else if(!strcmp(*p, "-t")){
            if(i+1 >= argc) break; 
            if(!strcmp(p[1], "max")){
                _param.max_trial_per_start = INT_MAX;
            }else{
                int ret = sscanf(p[1], "%d", &_param.max_trial_per_start);
                if(ret != 1) break;
            }
            i++, p++;
        }else if(!strcmp(*p, "-n")){
            if(i+1 >= argc) break;
             int ret = sscanf(p[1], "%d", &_param.ncores);
            if(ret != 1) break;
            i++, p++;
        }
        else if(!strcmp(*p, "-o")){
            if(i+1 >= argc) break;
            int ret = sscanf(p[1], "%d", &_param.n_violance_tolerance);
            if(ret != 1) break; 
            i++, p++;
        }else if(!strcmp(*p, "-m")){
            if(i+1 >= argc) break; 
            int ret = sscanf(p[1], "%f", &_param.multiplier);
            if(ret != 1) break;
            i++, p++;
        }else if(!strcmp(*p, "-v")){
            _param.verbose = 1;
        }else if(!strcmp(*p, "-maxsat")){
            _param.ismaxsat = 1;
        }else if(!strcmp(*p, "-wcnf")){
            _param.iswcnf = 1;
        }else if(!strcmp(*p,"-b")){
            if(i+1 >= argc) break;
            int ret = sscanf(p[1], "%f", &_param.beta);
            if(ret != 1) break;
            i++; p++;
        }else if(!strcmp(*p,"-a")){
            if(i+1 >= argc) break;
            int ret = sscanf(p[1], "%d", &_param.adapt);
            if(ret != 1) break;
            i++; p++;
        }  
        else if(i+1 == argc){
            _param.fin = fopen(*p, "r");
            if(!_param.fin){
                fprintf(stderr, "%s\n", strerror(errno));
                exit(1);
            }
            _param.fin_name = strdup(*p);
        }else{
            printf("Error: no such parameter\n");
            break;
        }
    }
    if(i != argc || !_param.fin){
        print_usage(argv[0], &_param);
        exit(0);
    }
    *param = _param;
}


int main(int argc, char **argv){
    Parameter param;
    get_parameter(argc, argv, &param);
    srand48(0);
    Formula formula(&param);
    formula.read_DIMACS();
    BDD bdd = BDD(&formula, param.L2);
    Optimizer_Portfolio op = Optimizer_Portfolio(formula.num_of_vars, &bdd, &param);
    op.solve();    
    return 0;
}

