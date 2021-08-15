import itertools
from boolean_formula import Formula
import sys
from pysat.card import *
from pysat.pb import *

def xor2cnf(xor):
    if len(xor)==3:
        return [[-xor[0],-xor[1],xor[2]], [xor[0],-xor[1],-xor[2]], [-xor[0],xor[1],-xor[2]],[xor[0],xor[1],xor[2]] ]
    if len(xor)==2:
        return [[xor[0],xor[1]],[-xor[0],-xor[1]]]
    
def de_xor(llist,start):
    oldstart = start
    exactFlag = 0
    if len(llist)<=3:
        return [list(llist),],start
    else:
        clist = []
        i = 0
        while True:
            if i>=len(llist):
                break
            elif i==len(llist)-1:
                exactFlag = 1
                break
            else:
                clist.append([-llist[i],llist[i+1],start+1])
                start+=1
                i+=2
    if exactFlag == 1:
        flist,_ = de_xor([llist[-1]] + list(range(oldstart+1,start+1)),start)
    else:
        flist,_ = de_xor(range(oldstart+1,start+1),start)
    return clist+flist,_

def kxor(n):
    return de_xor(range(1,n+1),n+1)

def xor100():
    f = open('benchmarks/sxor/x10.cnf','w+')
    n=10
    clauses,nv = kxor(n)
    f.write('p cnf '+repr(nv)+' '+repr(len(clauses))+'\n')
    for c in clauses:
        f.write('x ')
        for l in c:
            f.write(repr(l)+' ')
        f.write('0\n')

def sign(a):
    if a>0: return 1
    return -1

def pb_opb(filepath,topath):
    f = open(filepath,'r+')
    g = open(topath,'w+')
    flines = f.readlines()
    for line in flines:
        split = line.split()
        if len(split) == 0: continue
        if split[0] == 'p':
            n = int(split[2])
            m = int(split[3])
            g.write("* #variable= "+repr(n)+" #constraint= " +repr(m)+'\n')
        else:
            g.write(line)
    f.close()
    g.close()

    
def xor_blow(filepath,topath):
    encoding_no = 0
    f = open(filepath,'r+')
    g = open(topath,'w+')
    formula = Formula.read_DIMACS(filepath)
    new_formula = Formula()
    nv = len(formula.variables)
    nc = len(formula.clauses)
    ctype = formula._ctype
    start = nv
    clause_list = []
    m = 0
    for ci in range(len(formula.clauses)):
        if ctype[ci]=='c':
            if formula._klist[ci] == 1:
                clause_list.append([formula.clauses[ci]])
                m+=1
            else:
                cnf_from_pb = PBEnc.atleast(lits=formula.clauses[ci], weights=[1 for i in range(nv)],bound=abs(formula._klist[ci]),top_id = start, encoding=encoding_no)
                clause_list.append(cnf_from_pb.clauses)
                m += len(cnf_from_pb.clauses)
                if cnf_from_pb.nv > nv:
                    start = cnf_from_pb.nv
        if ctype[ci]=='p':
            cnf_from_pb = PBEnc.atleast(lits=formula.clauses[ci], weights=formula._coefs[ci] ,bound=abs(formula._klist[ci]),top_id = start, encoding=encoding_no)
            print('encoding cost '+repr(len(cnf_from_pb.clauses)))
            j = 0
            clause_list.append(cnf_from_pb.clauses)
            m += len(cnf_from_pb.clauses)
            if cnf_from_pb.nv > nv:
                start = cnf_from_pb.nv
        if ctype[ci]=='x':
            blow_clauses,_ = de_xor(formula.clauses[ci],start)
            start = _
            for b in blow_clauses:
                cnfs = xor2cnf(b)
                clause_list.append(cnfs)
                m += len(cnfs)
    g.write('p cnf '+repr(start)+' '+repr(m)+'\n')
    for cg in clause_list:
        for c in cg:
            for item in c:
                g.write(repr(item)+' ')
            g.write('0\n')
    f.close()
    g.close()

xor_blow(sys.argv[1],sys.argv[2])

