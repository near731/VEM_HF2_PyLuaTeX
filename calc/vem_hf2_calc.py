#Imports

import numpy as np
import sympy as sp
from sympy import init_printing
from sympy import *
from sympy import Matrix, symbols, solve_linear_system
from sympy.interactive import printing
from sympy.printing.latex import LatexPrinter, print_latex
from sympy import Eq
printing.init_printing(use_latex=True)

from vem_hf2_data import get_data
from mx_functions import ExtMatrix,SubMatrix,SubVector

#Code

code_2=2
code_3=2
code_4=2

def calculate():
    #get data

    data = get_data(code_2, code_3, code_4)

    a = data['a']
    m_0 = data['m_0']
    b = data['b']
    d = data['d']
    E = data['E']
    ρ = data['ρ']

    #2.feladat

    #Keresztmetszeti jellemzők 

    A=d**2*np.pi/4
    Iz=d**4*np.pi/64

    #Rudak hosszai

    L_1_a=b
    L_2_a=a

    #Merevségi mátrixok


    def K_matrix(L):
        K=Iz*E/(L**3)*sp.Matrix([[12,6*L,-12,6*L],
                              [6*L,4*L**2,-6*L,2*L**2],
                              [-12,-6*L,12,-6*L],
                              [6*L,2*L**2,-6*L,4*L**2]])
        return K

    #Elemenként

    K_1_a=K_matrix(L_1_a)
    K_2_a=K_matrix(L_2_a)
    

    #Globális merevségi mátrix

    #Elem-csomopont osszerendeles

    ecs_a = sp.Matrix([[1,2],[2,3]])

    eDOF1_a = [2*ecs_a[0,0]-1, 2*ecs_a[0,0], 2*ecs_a[0,1]-1, 2*ecs_a[0,1]]
    eDOF2_a = [2*ecs_a[1,0]-1, 2*ecs_a[1,0], 2*ecs_a[1,1]-1, 2*ecs_a[1,1]]

    K_glob_a = sp.Matrix(ExtMatrix(K_1_a, eDOF1_a, 6) + ExtMatrix(K_2_a, eDOF2_a, 6))
    
    #Tömeg mátrixok

    def M_i(L_i):
        M_i=ρ*A*L_i/420*sp.Matrix([[156,22*L_i,54,-13*L_i],
                           [22*L_i,4*L_i**2,13*L_i,-3*L_i**2],
                           [54,13*L_i,156,-22*L_i],
                           [-13*L_i,-3*L_i**2,-22*L_i,4*L_i**2]])
        return M_i

    #Elemként

    M_1_a=M_i(L_1_a)
    M_2_a=M_i(L_2_a)
    
    #Globális tömeg mátrix

    M_a = sp.Matrix(ExtMatrix(M_1_a, eDOF1_a, 6) + ExtMatrix(M_2_a, eDOF2_a, 6))

    #Kondenzált mennyiségek meghatározása

    #Kondenzált merevségi mátrix

    freeDoF_a=sp.Matrix([[3,4,6]])

    K_kond_a=SubMatrix(K_glob_a,freeDoF_a)

    #Kondenzált tömeg mátrix

    M_kond_a=SubMatrix(M_a,freeDoF_a)

    #Sajátkörfrekvenciák meghatározása:

    #Symbol

    α=sp.Symbol("α")

    #Alpha Calculation

    def alpha(K_i,M_i,n):
        A=K_i-α**2*M_i
        xarr = range(-n, n+1)
        yarr = [A.subs(α, x).det() for x in xarr]
        detA = expand(interpolating_poly(len(xarr), α, xarr, yarr))
        alphas=sp.solve(detA,α**2)
        α_vec = sp.Matrix([sp.re(α)**0.5 for α in alphas if sp.re(α) > 0 and abs(sp.im(α)) < 1e-5])
        α_vec = sorted(α_vec, key=lambda x: float(x))
        return α_vec

    #Sajátkörfrekvenciák

    α_a=alpha(K_kond_a,M_kond_a,3)
    α_1_a=α_a[0]
    α_2_a=α_a[1]
    α_3_a=α_a[2]

    #Sajátfrekveniák

    def freq(alpha_i):
        f_i=alpha_i/(2*np.pi)
        return f_i

    f_1_a=freq(α_1_a)
    f_2_a=freq(α_2_a)
    f_3_a=freq(α_3_a)

    #3.feladat

    #Rudak hosszai

    L_1_b=b/2
    L_2_b=b/2
    L_3_b=a

    #Merevségi mátrixok

    #Elemenként

    K_1_b=K_matrix(L_1_b)
    K_2_b=K_matrix(L_2_b)
    K_3_b=K_matrix(L_3_b)

    #Globális merevségi mátrix

    #Elem-csomopont osszerendeles

    ecs_b = sp.Matrix([[1,2],[2,3],[3,4]])

    eDOF1_b = [2*ecs_b[0,0]-1, 2*ecs_b[0,0], 2*ecs_b[0,1]-1, 2*ecs_b[0,1]]
    eDOF2_b = [2*ecs_b[1,0]-1, 2*ecs_b[1,0], 2*ecs_b[1,1]-1, 2*ecs_b[1,1]]
    eDOF3_b = [2*ecs_b[2,0]-1, 2*ecs_b[2,0], 2*ecs_b[2,1]-1, 2*ecs_b[2,1]]

    K_glob_b = sp.Matrix(ExtMatrix(K_1_b, eDOF1_b, 8) + ExtMatrix(K_2_b, eDOF2_b, 8)+ ExtMatrix(K_3_b, eDOF3_b, 8))
    

