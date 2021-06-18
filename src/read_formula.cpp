#include "../include/read_formula.h"
#include<string>
#include<fstream>
#include<iostream>
#include<sstream>
#include<stdlib.h>
#include <iterator>
#include <math.h>

Formula::Formula(Parameter *param){
    this->num_of_vars = 0;
    this->clauses = new std::vector<std::vector<int> *>;
    this->coefs = new std::vector<std::vector<int> *>;
    this->clause_weights = new std::vector<double>;
    this->clause_weights_original = new std::vector<double>;
    this->clause_type = new std::vector<char>;
    this->klist = new std::vector<int>;
    this->comparators = new std::vector<int>;
    this->param = param;
}

void Formula::add_clause(std::vector<int>* literals, int k, char ctype, double weight, std::vector<int> *coefs = NULL, int comparator =0){
   this->clauses->push_back(literals);
   this->klist->push_back(k);
   this->clause_type->push_back(ctype);
   this->clause_weights->push_back(weight);
   this->clause_weights_original->push_back(weight);
   this->coefs->push_back(coefs);
   this->comparators->push_back(comparator);
}

void Formula::print(){
    for (int i=0;i<this->clauses->size();i++){
        std::cout<<(*this->clause_type)[i]<<" ";
        for (int j=0;j<(*this->clauses)[i]->size();j++){
            std::cout<<(*(*this->clauses)[i])[j]<<" ";
        }
        std::cout<<std::endl;
    }
}

static void PB_canonicalize(std::vector<int> *literals, std::vector<int> *coefs, int *k, int *comparator){ // comparator: >= (1), =(2), <= (3)
    if (*comparator == 3){
        *comparator = 1;
        *k = -(*k);
        for ( int i = 0; i < coefs->size(); i++){
            (*coefs)[i] = -((*coefs)[i]);
        }
    }
    for ( int i = 0; i < coefs->size(); i++){
        if ((*coefs)[i] < 0 ){
            (*coefs)[i] = -((*coefs)[i]);
            (*literals)[i] = -(*literals)[i];
            (*k) += (*coefs)[i];
        }
    }       
}

void Formula::read_DIMACS(){
    std::string file = this->param->fin_name;
    std::ifstream input_file;
    input_file.open(file.c_str());
    if(!input_file.is_open()){
        std::perror("error open");
        exit(1);
    }
    std::string line;
    while(std::getline(input_file,line)){
        std::stringstream ss(line);
    	std::istream_iterator<std::string> begin(ss);
    	std::istream_iterator<std::string> end;
    	std::vector<std::string> split(begin, end);
        std::vector<int> *literals = new std::vector<int>;
        int weight;
        char ctype;
        if ((split.size() >= 2) && (split[1] == "#variable=") ){
            this->num_of_vars = std::stoi(split[2]);
        }
        if ( (split.size()==0) || (split[0]=="c") || (split[0][0])=='*' || (split[0][0]=='c')) continue;
        if ( split[0] == "p"){
            this->num_of_vars = std::stoi(split[2]);   
            if (split[1] == "wcnf"){
                this->param->iswcnf = 1;
                this->param->wcnf_weight = std::stoi(split[4]);
            }  
        }
        else if (split[0] == "g"){
            std::vector<int> *coefs = new std::vector<int>(this->num_of_vars,1);
            int k = std::stoi(split[1]);
            weight = this->compute_clause_weight(this->num_of_vars,k,'g');
            if (k>=0){
                for(int i=1; i<=this->num_of_vars;i++){
                    literals->push_back(i);}
                this->add_clause(literals,k,'c',weight, coefs);
	    }
            else{
                for(int i=1; i<=this->num_of_vars;i++){
                    literals->push_back(-i);}
                this->add_clause(literals, this->num_of_vars + k, 'c', weight, coefs);
            }
        }
        else if (split[0] == "x"){
            for(int i=1; i<split.size()-1;i++){
                literals->push_back(std::stoi(split[i]));
            }
            weight = this->compute_clause_weight(literals->size(),1,'x');
            this->add_clause(literals, 1 , 'x', weight);
        }
        else if (split[0] == "d" ){
            int k = std::stoi(split[1]);
            if (k>=0){
            	for(int i=2; i<split.size()-1;i++){
                	literals->push_back(std::stoi(split[i]));
            	}
                weight = this->compute_clause_weight(literals->size(),k,'d');
                std::vector<int> *coefs = new std::vector<int>(literals->size(),1);
                this->add_clause(literals, k , 'c', weight,coefs);
            }
             else{
                for(int i=2; i<split.size()-1;i++){
                        literals->push_back(-std::stoi(split[i]));
                }
                weight = this->compute_clause_weight(literals->size(),k,'d');
                std::vector<int> *coefs = new std::vector<int>(literals->size(),1);
                this->add_clause(literals, literals->size() + k , 'c', weight, coefs);
            }
        }
        else if (split[0] == "n" ){
           for(int i=1; i<split.size()-1;i++){
                literals->push_back(std::stoi(split[i]));
            }
            weight = this->compute_clause_weight(literals->size(),1,'n');
            this->add_clause(literals, 1 , 'n', weight);            
        }
        else if (split[1][0] == 'x' ){  // pseudo-Boolean constraints
            std::vector<int> *coefs = new std::vector<int>;
            for(int i=0; i < (split.size() - 3) / 2;i++){
                coefs->push_back(std::stoi(split[i*2]));
                literals->push_back(std::stoi(split[i*2+1].erase(0,1)));
            }
            std::string comparator_string = split[split.size()-3];
            int comparator = 0;
            int k = std::stoi(split[split.size() -2]);
            if(comparator_string == ">=") comparator = 1;
            else if(comparator_string == "=") comparator = 2;
            else if(comparator_string == "<=") comparator = 3;
            PB_canonicalize(literals, coefs, &k, &comparator); // comparator: >= (1), =(2) 
            weight = this->compute_clause_weight(literals->size(),k,'p');
            this->add_clause(literals, k , 'p', weight,coefs,comparator);
        }
        else{
            if ( this->param->ismaxsat && this->param->iswcnf){
                  // read wcnf
                  weight = std::stoi(split[0]);
                  for(int i=1; i<split.size()-1;i++){
                    literals->push_back(std::stoi(split[i]));
                 }
            }
            else{
                for(int i=0; i<split.size()-1;i++){
                    literals->push_back(std::stoi(split[i]));
    	        }
                weight = this->compute_clause_weight(literals->size(),1,'c');
            }
            std::vector<int> *coefs = new std::vector<int>(literals->size(),1);
            this->add_clause(literals, 1 , 'c', weight, coefs);
        }
    }
}


inline int min(int a,int b){
    return (a < b)?a:b;
}

double Formula::compute_clause_weight(int n, int k, char ctype){
    return n;
}
