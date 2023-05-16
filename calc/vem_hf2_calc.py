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

from .vem_hf2_data import get_data
from .mx_functions import ExtMatrix,SubMatrix,SubVector

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

    #Tömeg mátrixok

    #Elemként

    M_1_b=M_i(L_1_b)
    M_2_b=M_i(L_2_b)
    M_3_b=M_i(L_3_b)

    #Globális tömeg mátrix

    M_b = sp.Matrix(ExtMatrix(M_1_b, eDOF1_b, 8) + ExtMatrix(M_2_b, eDOF2_b, 8)+ ExtMatrix(M_3_b, eDOF3_b, 8))

    #Kondenzált mennyiségek meghatározása

    #Kondenzált merevségi mátrix

    freeDoF_b=sp.Matrix([[3,4,5,6,8]])

    K_kond_b=SubMatrix(K_glob_b,freeDoF_b)

    #Kondenzált tömeg mátrix

    M_kond_b=SubMatrix(M_b,freeDoF_b)

    #Sajátkörfrekvenciák meghatározása:

    #Sajátkörfrekvenciák

    α_b=alpha(K_kond_b,M_kond_b,4)
    α_1_b=α_b[0]
    α_2_b=α_b[1]
    α_3_b=α_b[2]

    #Sajátfrekvenciák

    f_1_b=freq(α_1_b)
    f_2_b=freq(α_2_b)
    f_3_b=freq(α_3_b)

    #4.feladat

    #Rudak hosszai

    L_1_c=L_1_b
    L_2_c=L_2_b
    L_3_c=L_3_b

    #Merevségi mátrixok

    K_glob_c=K_glob_b

    #Tömeg mátrixok

    M_k=sp.Matrix([[0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0],
                [0,0,0,0,m_0,0,0,0],
                [0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0]])

    M_c=M_b+M_k

    #Mennyiségek kondenzálása

    freeDoF_c=freeDoF_b

    #Kondenzált merevségi mátrix

    K_kond_c=K_kond_b

    #Kondenzált tömeg mátrix

    M_kond_c=SubMatrix(M_c,freeDoF_c)

    #Sajátkörfrekvenciák meghatározása:

    #Sajátkörfrekvenciák

    α_c=alpha(K_kond_c,M_kond_c,5)

    α_1_c=α_c[0]
    α_2_c=α_c[1]
    α_3_c=α_c[2]

    #Sajátfrekvenciák

    f_1_c=freq(α_1_c)
    f_2_c=freq(α_2_c)
    f_3_c=freq(α_3_c)

    #5.feladat

    #Rudak hosszai

    L_1_d=L_1_b
    L_2_d=L_2_b
    L_3_d=L_3_b

    #Merevségi mátrixok

    K_glob_d=K_glob_b

    #Tömeg mátrixok

    def M_e(L_i):
        M_e=ρ*A*L_i/2*sp.Matrix([[1,0,0,0],
                             [0,L_i**2/12,0,0],
                             [0,0,1,0],
                             [0,0,0,L_i**2/12]])
        return M_e

    #eDOFs

    eDOF1_d=eDOF1_b
    eDOF2_d=eDOF2_b
    eDOF3_d=eDOF3_b

    #Tömeg mátrixok

    M_1_d=M_e(L_1_d)
    M_2_d=M_e(L_2_d)
    M_3_d=M_e(L_3_d)

    M_d = sp.Matrix(ExtMatrix(M_1_d, eDOF1_d, 8) + ExtMatrix(M_2_d, eDOF2_d, 8)+ ExtMatrix(M_3_d, eDOF3_d, 8))
    M_d+=M_k

    #Mennyiségek kondenzálása

    freeDoF_d=freeDoF_b

    #Kondenzált merevségi mátrix

    K_kond_d=K_kond_b

    #Kondenzált tömeg mátrix

    M_kond_d=SubMatrix(M_d,freeDoF_d)

    #Sajátkörfrekvenciák meghatározása:

    #Sajátkörfrekvenciák

    α_d=alpha(K_kond_d,M_kond_d,5)
    α_1_d=α_d[0]
    α_2_d=α_d[1]
    α_3_d=α_d[2]

    #Sajátfrekvenciák

    f_1_d=freq(α_1_d)
    f_2_d=freq(α_2_d)
    f_3_d=freq(α_3_d)

    V={
       "A":A,
       "Iz":Iz,

       "L_1_a":L_1_a,
       "L_2_a":L_2_a,

       "K_1_a":K_1_a,
       "K_2_a":K_2_a,
       "K_glob_a":K_glob_a,

       "M_1_a":M_1_a,
       "M_2_a":M_2_a,
       "M_a":M_a,

       "K_kond_a":K_kond_a,
       "M_kond_a":M_kond_a,

        "α_1_a":α_1_a,
        "α_2_a":α_2_a,
        "α_3_a":α_3_a,

        "f_1_a":f_1_a,
        "f_2_a":f_2_a,
        "f_3_a":f_3_a,

        "L_1_b":L_1_b,
        "L_2_b":L_2_b,
        "L_3_b":L_3_b,

       "K_1_b":K_1_b,
       "K_2_b":K_2_b,
       "K_3_b":K_3_b,
       "K_glob_b":K_glob_b,

       "M_1_b":M_1_b,
       "M_2_b":M_2_b,
       "M_3_b":M_3_b,
       "M_b":M_b,

       "K_kond_b":K_kond_b,
       "M_kond_b":M_kond_b,

        "α_1_b":α_1_b,
        "α_2_b":α_2_b,
        "α_3_b":α_3_b,

        "f_1_b":f_1_b,
        "f_2_b":f_2_b,
        "f_3_b":f_3_b,

        "M_c":M_c,
        "M_kond_c":M_kond_c,

        "α_1_c":α_1_c,
        "α_2_c":α_2_c,
        "α_3_c":α_3_c,

        "f_1_c":f_1_c,
        "f_2_c":f_2_c,
        "f_3_c":f_3_c,

        "M_1_d":M_1_d,
        "M_2_d":M_2_d,
        "M_3_d":M_3_d,
        "M_d":M_d,

        "M_kond_d":M_kond_d,

        "α_1_d":α_1_d,
        "α_2_d":α_2_d,
        "α_3_d":α_3_d,

        "f_1_d":f_1_d,
        "f_2_d":f_2_d,
        "f_3_d":f_3_d,
    }
    for k, v in V.items():
        data[k] = v
    return data