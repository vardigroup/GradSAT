#include <dlib/optimization.h>
#include <iostream>
#include "../include/optimize.h"
#include <random>
#include <algorithm>
#include <stdlib.h>
#define CHANGE_WEIGHTS
#include "omp.h"
using namespace std;
using namespace dlib;

static double sum(std::vector<double> *v){
    double s = 0;
    for (int i = 0; i < v->size(); i++){
        s += (*v)[i];
    }
    return s;
}

Optimizer::Optimizer(){
}

Optimizer::Optimizer(int n, BDD *bdd, std::string name, Parameter *param){
    srand(time(NULL));
    this->bdd = bdd;
    this->num_of_vars = n;
    this->solved_flag = 0;
    this->x = new std::vector<double>;
    this->unsat_clauses = new std::vector<int>;
    for(int i=0; i< n; i++){
        double xi = (double) std::rand();
        this->x->push_back( 1 - 2 * (xi/(double)RAND_MAX));
    }
    this->optimizer_name = name;
    this->param = param;
    this->number_of_trials_this_start = 0;
}

double Optimizer::fval_for_optimizer(const column_vector& m)
{
    std::vector<double> x;
    for (int i=0; i < this->bdd->num_of_vars; i++){
        x.push_back((double)m(i));
    }
    return this->bdd->fval(&x);
}

const column_vector Optimizer::grad_for_optimizer (const column_vector& m)
{
    std::vector<double> *grad;
    std::vector<double> x;
    for (int i=0; i < this->bdd->num_of_vars; i++){
        x.push_back((double)m(i));
    }
    column_vector res(this->bdd->num_of_vars);
    grad = this->bdd->grad(&x);
    for (int i=0; i < this->bdd->num_of_vars; i++){
        res(i) = (*grad)[i];
    }
    return res;
}

double fval_for_optimizer_nlopt(const std::vector<double> &x, std::vector<double> &grad, void *data){
    Optimizer *opt = (Optimizer *) data;
    std::vector<double> *grad_p = opt->bdd->grad(&x);
    for ( int i = 0; i<grad.size(); i++) grad[i] = (*grad_p)[i];
    return opt->bdd->fval(&x);
}

double fval_for_optimizer_nlopt_gdfree(const std::vector<double> &x, std::vector<double> &grad, void *data){
    Optimizer *opt = (Optimizer *) data;
    return opt->bdd->fval(&x);
}

void Optimizer::update_weights(){
    for ( int i = 0; i < this->unsat_clauses->size(); i++){
        (*this->bdd->clause_weights) [ (*this->unsat_clauses)[i] ] *= this->param->multiplier;
    }
    this->bdd->sum_of_clause_weights = sum(this->bdd->clause_weights);
}

static std::vector<double> *rounding( column_vector *m, int n){
    std::vector<double> *res = new std::vector<double>;
    for ( int i = 0; i < n; i++){
        if ( (*m)(i) > 0 )
            res->push_back(1);
        else res->push_back(-1);
    }
    return res;  
}

static std::vector<double> *rounding( std::vector<double> *x, int n){
    std::vector<double> *res = new std::vector<double>;
    for ( int i = 0; i < n; i++){
        if ( (*x)[i] > 0 )
            res->push_back(1);
        else res->push_back(-1);
    }
    return res;
}

static void random_restart(std::vector<double> *x){
    int n = x->size();
    for(int i=0; i< n; i++){
        double xi = (double) std::rand();
        (*x)[i] =  1 - 2 * (xi/(double)RAND_MAX);
    }

}

static int count(std::vector<double> *x){
    int count = 0;
    for( int i =0; i< x->size();  i++){
        if ((*x)[i] < 0) count++;
    }
    return count;
}

static int count(column_vector &m, int n){
    int count = 0;
    for( int i =0; i< n;  i++){
        if (m(i) < 0) count++;
    }
    return count;
}

static void print_vector(std::vector<double> *x){
    int n = x->size();
    for(int i=0; i< n; i++){
        std::cout<<(*x)[i]<<",";
    }
    std::cout<<std::endl;
}

static void print_solution(std::vector<double> *x){
    int n = x->size();
    std::cout<<"v ";
    for(int i=0; i< n; i++){
        if ((*x)[i] < 0){
            std::cout<<i+1<<" ";
        }
        else 
            std::cout<<-(i+1)<<" ";
    }
    std::cout<<std::endl;
}

static void print_solution(column_vector &m, int n){
    std::cout<<"v ";
    for(int i=0; i< n; i++){
        if (m(i) < 0){
            std::cout<<i+1<<" ";
        }
        else
            std::cout<<-(i+1)<<" ";
    }
    std::cout<<std::endl;
}

double Optimizer::minimize()
{
    random_restart(this->x);
    int n = this->num_of_vars;
    static int min_unsat_for_maxsat = 1e9;
    double tol = 1e-40;
    while(1){
        if (  number_of_trials_this_start  >= this->param->max_trial_per_start ){
          //   if (this->param->verbose) std::cout<<"c Random restart"<<std::endl;
             random_restart(this->x);
             this->bdd->restart_weights();
             this->number_of_trials_this_start = 0;
        }
        this->number_of_trials_this_start ++;
        double stop_val = this->bdd->sum_of_clause_weights;
        column_vector starting_point(n);
        std::vector<double> *x = new std::vector<double>(n);
        for(int i = 0; i < n; i++){
            (*x)[i] = (*this->x)[i];
            starting_point(i) = (double) (*x)[i];
        }
        std::vector<double> *roundedx;
        std::string package;
        if (this->optimizer_name=="CG"){
            if (this->param->verbose) 
                std::cout<<"c Solver: CG"<<std::endl;
            package = "DLIB";
            find_min_box_constrained(cg_search_strategy(), 
                            objective_delta_stop_strategy(tol),  
                           [this](const column_vector& a){ return this->fval_for_optimizer(a);}, 
                           [this](const column_vector& a){ return this->grad_for_optimizer(a);},
                           starting_point, -1.0, 1.0);
            roundedx = rounding(&starting_point,n);
                std::cout<<"c CG finds local maximum with value "<<-this->bdd->fval(roundedx)
                    <<". Goal: "<<this->bdd->sum_of_clause_weights<< std::endl;
        }
        else if (this->optimizer_name=="BFGS"){
            if (this->param->verbose) 
                std::cout<<"c Solver: BFGS"<<std::endl;
            package = "DLIB";
            find_min_box_constrained(bfgs_search_strategy(),
                            objective_delta_stop_strategy(tol),
                           [this](const column_vector& a){ return this->fval_for_optimizer(a);},
                           [this](const column_vector& a){ return this->grad_for_optimizer(a);},
                           starting_point, -1.0, 1.0);
            roundedx = rounding(&starting_point,n);
            std::cout<<"c BFGS finds local maximum with value "<<-this->bdd->fval(roundedx)
                    <<". Goal: "<<this->bdd->sum_of_clause_weights<< std::endl;
        }
        else if ( this->optimizer_name=="SLSQP"){
            if (this->param->verbose) 
                std::cout<<"c Solver: SLSQP"<<std::endl;
            package = "NLOPT";
            nlopt::opt opt(nlopt::LD_SLSQP, n);
            std::vector<double> lb(n);
            std::vector<double> ub(n);
            for ( int i =0; i<n;i++){
                lb[i] = -1;
                ub[i] = 1;
            }
            opt.set_lower_bounds(lb);
            opt.set_upper_bounds(ub);
            opt.set_min_objective(fval_for_optimizer_nlopt, this);
            opt.set_xtol_rel(tol);
            opt.set_ftol_abs(tol);
            double minf;
            try{
                nlopt::result result = opt.optimize(*x, minf);
                    std::cout << "c SLSQP founds local maximum with value"<< std::setprecision(10) << -minf 
                    <<". Goal: "<<this->bdd->sum_of_clause_weights<< std::endl;
                }
            catch(std::exception &e) {
                if (this->param->verbose) 
                    std::cout << "c nlopt failed: " << e.what() << std::endl;
            }
            roundedx = rounding(x,n);
        }
        else if ( this->optimizer_name=="MMA"){
            if (this->param->verbose) 
                std::cout<<"c Solver: MMA"<<std::endl;
            package = "NLOPT";
       // nlopt::opt opt(nlopt::LD_MMA, n);
            nlopt::opt opt(nlopt::LD_CCSAQ, n);
            std::vector<double> lb(n);
            std::vector<double> ub(n);
            for ( int i =0; i<n;i++){
                lb[i] = -1;
                ub[i] = 1;
            }
            opt.set_lower_bounds(lb);
            opt.set_upper_bounds(ub);
            opt.set_min_objective(fval_for_optimizer_nlopt, this);
            opt.set_xtol_rel(tol);
            opt.set_ftol_abs(tol);
            double minf;
            try{
                nlopt::result result = opt.optimize(*x, minf);
                    std::cout << "c MMA finds local minimum with value "
                    << std::setprecision(10) << -minf 
                    <<". Goal: "<<this->bdd->sum_of_clause_weights<< std::endl;
            }
            catch(std::exception &e) {
                if (this->param->verbose) 
                    std::cout << "c nlopt failed: " << e.what() << std::endl;
            }
            roundedx = rounding(x,n);
        }
        else{ std::cout<<"c undefined optimizer: "<<this->optimizer_name<<std::endl; exit(0);}
        
        
        int num_unsat_clause, unsat_weight;
        this->bdd->verify_solution(roundedx,this->unsat_clauses, &num_unsat_clause, &unsat_weight);
        //double fval = fval_for_optimizer(starting_point);
        #ifdef CHANGE_WEIGHTS
            this->update_weights();
        #endif
        if (this->param->ismaxsat && (!this->param->iswcnf)){
            #pragma omp critical
            {
            if (num_unsat_clause < min_unsat_for_maxsat){
                min_unsat_for_maxsat = num_unsat_clause;
             }
                std::cout<<"o "<<min_unsat_for_maxsat<<std::endl;
            }
        }
        else if (this->param->ismaxsat && (this->param->iswcnf)){
            #pragma omp critical
            {
            if ( (unsat_weight < this->param->wcnf_weight) && (num_unsat_clause < min_unsat_for_maxsat)){
                min_unsat_for_maxsat = num_unsat_clause;
             }
                std::cout<<"o "<<min_unsat_for_maxsat<<std::endl;
            }
        }
        else{
           #pragma omp critical
            {
            if (num_unsat_clause < min_unsat_for_maxsat){
                std::cout<<"o "<<num_unsat_clause<<std::endl;
                min_unsat_for_maxsat = num_unsat_clause;
             }
            }
        }
        #pragma omp critical
        {
        if (num_unsat_clause<=this->param->n_violance_tolerance){ 
            cout<<"c solved by "<<this->optimizer_name<<". trials: "<<this->number_of_trials_this_start<<std::endl;
            cout<<"s SATISFIABLE "<<std::endl;
            if ( package == "DLIB" ){ 
                    print_solution(starting_point,n);
            }
            else if ( package == "NLOPT"){
                    print_solution(x);
            }
            this->solved_flag = 1;
            exit(0);
      }  
      }
    }
    return 0;
}

Optimizer_Portfolio::Optimizer_Portfolio(){
}


Optimizer_Portfolio::Optimizer_Portfolio(int num_of_vars, BDD *original_bdd, Parameter *param){
    this->num_of_vars = num_of_vars;
    this->bdd = original_bdd;
    this->param = param;
}

void Optimizer_Portfolio::solve(){
    BDD bdd_group[this->param->ncores];
    Optimizer optimizer_group[this->param->ncores];
    std::vector<std::string> solvers{"CG","SLSQP","MMA","BFGS"};
    if (this->param->solver != "PORTFOLIO"){
       solvers.clear();
       solvers.push_back(this->param->solver);
       cout<<"solver: "<<this->param->solver<<endl;
    }
    for ( int i = 0; i < this->param->ncores; i++){
        bdd_group[i] = BDD(this->bdd);
        optimizer_group[i] = Optimizer(this->num_of_vars, &bdd_group[i], solvers[i % solvers.size()], this->param);
    }
    int num_trials = 0;
    bool solved_flag = 0;
    if (this->param->verbose)
        cout<<"c ncores "<<this->param->ncores<<endl;
    int nthreads =  param->ncores;
    #pragma omp parallel for num_threads(this->param->ncores)
    for ( int i = 0; i < nthreads; i++){
            int tid = omp_get_thread_num();
            optimizer_group[i].minimize();
            solved_flag |= optimizer_group[i].solved_flag;                
    }
    cout<<"c number of trials: "<<num_trials<<endl;
}
