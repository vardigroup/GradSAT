#ifndef MAIN_H
#define MAIN_H
#include <random>
#include <algorithm>
struct Parameter {
    std::string solver;
    int k;
    double eps;
    int max_iter;
    double multiplier;
    int verbose;
    FILE *fin;
    char *fin_name;
    double beta;
    int adapt;
    int ncores;
    int max_trial_per_start;
    float L2;
    int n_violance_tolerance;
    int ismaxsat;
    int iswcnf;
    int wcnf_weight;
};

#endif
