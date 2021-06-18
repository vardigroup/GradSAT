#ifndef BDD_GRADIENT_H
#define BDD_GRADIENT_H
#include "../include/read_formula.h"
#include<string>
#include "cudd.h"
#include <sys/types.h>
#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <map>
#include <iostream>


class BDD{
    public:
    DdManager *gbm;
    std::vector<double> *clause_weights;
    std::map<DdNode*,double> *forward_message;
    std::map<DdNode*, double> *backward_message;  
    std::map<DdNode*, std::vector<DdNode*> * > *hipar;  
    std::map<DdNode*, std::vector<DdNode*> * > *lopar;  
    std::vector<DdNode*> *roots;
    int num_of_vars;
    double sum_of_clause_weights;
    Formula *formula;
    float L2;    
    
    BDD();
    BDD(Formula *formula, float L2);
    BDD(BDD *bdd);
    void change_weights_in_BDD(std::vector<double> *clause_weights);
    void restart_weights();
    std::vector<double> *grad(const std::vector<double> *x);
    double fval(const std::vector<double> *x);
    void print_dd (DdManager *gbm, DdNode *dd, int n, int pr );   
    void verify_solution(std::vector<double> *x, std::vector<int> *unsat_clauses, int *num_unsat_clauses, int *unsat_weight);
    
    private:
    int generate_parents;
    void message_clean();
    void backward_message_clean();
    void build_BDD_for_clause(Formula *formula, int ci);
    void forward_pass(std::vector<double> *x);
    void backward_pass(std::vector<double> *x, std::vector<double> *grad);
    void combine_clauses(int k);
    DdNode *BDD_for_PB(DdManager *gbm, std::vector<int> *literals, std::vector<int> *coefs, int rhs, int size, int sum, int material_left, int comparator, std::map<std::pair<int,int>, DdNode *> *nodes_storage);
    DdNode *BDD_for_XOR(DdManager *gbm, std::vector<int> *literals, int size, int product, std::map<std::pair<int,int>, DdNode *> *nodes_storage);
};
#endif
