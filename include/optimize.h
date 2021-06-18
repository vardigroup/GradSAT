#ifndef OPTIMIZE_H
#define OPTIMIZE_H
#include<vector>
#include<string>
#include<dlib/optimization.h>
#include<iostream>
#include "bdd_gradient.h"
#include <math.h>
#include <nlopt.hpp>
#include <iomanip>
#include "main.h"
typedef dlib::matrix<double,0,1> column_vector;

class Optimizer{
    public:
    int num_of_vars;
    BDD *bdd; 
    std::vector<double> *x; 
    void update_weights();    
    std::vector<int> *unsat_clauses;   
    Optimizer();
    Optimizer(int n, BDD *bdd,std::string name, Parameter *param);
    double minimize();
    double fval(const column_vector& m);
    Parameter *param;
    column_vector a_grad (const column_vector& m);
    bool solved_flag;
    int verbose;    
        int number_of_trials_this_start; 
        std::string optimizer_name;
    private:
    double fval_for_optimizer(const column_vector& m);
    const column_vector grad_for_optimizer (const column_vector& m);
};


class Optimizer_Portfolio{
    public:
        int num_of_vars;
        std::string optimizer_name;
        Parameter *param;
        BDD *bdd;
        Optimizer_Portfolio();
        Optimizer_Portfolio(int num_of_vars, BDD *original_bdd, Parameter *param);
        void solve();
};
#endif
