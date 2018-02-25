#libraries import
from numpy import sqrt, sin, cos, pi
from pylab import plot, show, bar
import numpy as np
import matplotlib.pyplot as plt
from scipy import *
import scipy.special
from numpy.linalg import *
from sympy import *
import sympy as sp
import sympy as sp

t = sp.symbols('t')

def plotter(array,l,u,n,t):
    size = (u-l)/n
    t = float(t)
    pl = [0 for i in range(0, n+1)]
    for row in range(1,n+1):
        pl[row] = l+size*row
        if(t> pl[row-1] and t<=pl[row]):
            print(array[row])
    plt.step(pl, array)
    plt.show()

def answeri(array,l,u,n,t):
    size = (u - l) / n
    t = float(t)
    pl = [0 for i in range(0, n)]
    for row in range(0, n):
        pl[row] = l + size * row
    for row in range(0, n-1):
        if (t >= pl[row] and t < pl[row+1]):
            return array[row]
    if(t>=pl[n-1]):
        return array[n-1]
#defining haar
def haar(t):
     return sp.Piecewise((1,(t>=0)&(t<0.5)),(-1,(t>=0.5)&(t<1)),(0,(t>=1)|(t<0)))


#defining h[i]s
h=[0 for i in range(0,10000)]
h[0]=1
#k is differnt  from k_in
for j in range(0, 5):
    for k in range(0, (2 ** j) - 1 + 1):
        i = (2 ** j) + k
        h[i] = (2 ** (j / 2)) * haar((2 ** j) * t - k)

def H_k_particular(t_value, k):
    matrix = [0 for x in range(k)]
    for row in range(k):
        matrix[row] = float(h[i].subs(t,t_value))
    return matrix

def H_k_general(k):
    matrix = [0 for x in range(k)]
    for row in range(k):
        matrix[row] = h[i]
    return matrix


#a function which returns Hk*k
def Hkk_m(k):
    Hkk = [[0 for x in range(k)] for y in range(k)]
    for i in range (0,k):

        for j in range (0,k):
            if(i==0):
                Hkk[i][j] = 1
            else:
                Hkk[i][j] = float(h[i].subs(t,(j+1/2)/k))
    return Hkk

#Epsilon function required to create F_alpha
def Epsilon(k):
    return (k+1) ** (alpha_input + 1) - 2 * (k ** (alpha_input + 1)) + (k - 1) ** (alpha_input + 1)

#function as given in the example
def func_A(t):
    val =t**(1/3)
    return val

#function as given in the example
def Func_Example(t):
    return 3/(2*scipy.special.gamma(2/3))*(t**(2/3)) + t**(4/3)

#define y=x
def func_u(tt):
    return tt

#function to make recursive QKK for calculating error
def recursive_matrix(k):
    if(k==1):
        m = np.ones((1,1))
        m[0][0]=(1/2)
        return m
    matrix = [[0 for x in range(k)] for y in range(k)]
    qk2k2 = recursive_matrix(floor(k / 2))
    for i in range (0,floor(k/2)):
        for j in range (0,floor(k/2)):
            matrix[i][j]  = qk2k2[i][j]
    hk2k2 = Hkk_m(floor(k / 2))
    for i in range (floor(k/2),k):
        for j in range (0,floor(k/2)):

            matrix[i][j]  = -1*hk2k2[i-floor(k/2)][j]
    hk2k2in = inv(Hkk_m(floor(k / 2)))
    for i in range(0,floor(k/2)):
        for j in range(floor(k/2), k ):

            matrix[i][j] =  hk2k2in[i][j-floor(k/2)]
    return matrix

#taking input k
k_input = int(input("Enter value of k:"))
alpha_input = float(input("Enter value of alpha:"))

#initializing Hkk based on k_input
Hkk = Hkk_m(k_input)

#initializing _alpa_intermediate based on k_input and alpha
#f_alpha intermediate is the matrix consisting of epsilon functions 0s and 1s
F_Alpha_inter = np.zeros((k_input, k_input))
for row in range(k_input):
    counter=1
    for col in range (row,k_input):
        if col==row:
            F_Alpha_inter[row][col]=1
        else:
            F_Alpha_inter[row][col]=Epsilon(counter)
            counter = counter + 1

#finding the final F_Alpha
F_Alpha = 1 / ((k_input ** alpha_input) * (scipy.special.gamma(alpha_input + 2))) * F_Alpha_inter
Hkk_inverse = inv(Hkk)

P_Alpha = np.dot(Hkk, F_Alpha)
P_Alpha = np.dot(P_Alpha, Hkk_inverse)

#initializing the matrix f
f_k=np.zeros((k_input,1))
for i in range(k_input):
    f_k[i]=Func_Example((i + 0.5) / k_input)

A_kk=np.zeros((k_input,k_input))
for row in range(k_input):
     for col in range(k_input):
         if col==row:
              A_kk[row][col]=func_A((row+0.5)/k_input)

#initializing the c matrix used as cT
c_k=np.ones((k_input,1))


# intermediate value to calculate ctranspose
Finder= np.ones((k_input, k_input))
#temp2 is H
# temp2=[[0 for x in range(k_input)] for y in range(k_input)]
# for row in range(0, k_input):
#     for col in range(0, k_input):
#         if (row == 0):
#             Hkk[row][col] = 1
#         else:
#             Hkk[row][col] = float(h[row].subs(t, (col + 0.5) / k_input))
#intermediate value to calculate cT
inter_val=np.dot(P_Alpha, Hkk)
inter_val = np.dot(inter_val,A_kk)

Finder = inter_val + Hkk
Finder_transpose = np.transpose(Finder)

#linalg.solve(A,B) -- AX=B -- returns X
c_k = np.linalg.solve(Finder_transpose, f_k)
c_k_transpose=np.transpose(c_k)

Hkt = H_k_general(k_input)
c_k = np.dot(c_k_transpose,P_Alpha)
ut_general = np.dot(c_k, Hkk)
ut_general_func = ut_general[0]

utccollect =np.insert(ut_general[0],0,0.0)
t_input = float (input("Enter the value of t"))
#plotting the function a well as checking the estimated value
plotter(utccollect,0,1,k_input,t_input)

# h_k_col=np.zeros((k,1))
# for i in range(k):
#     h_k_col[i]=ttt[i]
#
# c_k_transpose=np.transpose(c_k)
# jjj=np.dot(temp,h_k_col)
# ans= np.dot(c_k_transpose, jjj)
# kkk=ans[0]
# ans1=kkk[0]
# err=math.fabs(ans1-t_in)
# print  (ans1, " Delta ERROR  :  ",err)

QKK = recursive_matrix(k_input)

m = 1

#Qkkm to calculate Q to the power m
Qkkm = [[0 for x in range(k_input)] for y in range(k_input)]
Qkkm = QKK

for j in range(1,m):
    Qkkm = np.dot(Qkkm,QKK)

#this is the transpose of ut
ut_tr = [0 for x in range(k_input)]
#ut is the matrix ut of error
#ut_general func is the calculated estimated function
for row in range (0,k_input):
    ut_tr[row] = answeri(ut_general_func,0,1,k_input,(row+1/2)/k_input)
ut = np.transpose(ut_tr)


bi =[0 for x in range(k_input)]
#to calculate bis of the equation utranspose = btranspose * qkk * hkk
Dot_of_qkk_hkk = np.dot(Qkkm, Hkk)
dotofqkkhkk_transpose = np.transpose(Dot_of_qkk_hkk)

bi = np.linalg.solve(dotofqkkhkk_transpose,ut)
bi_transpose = np.transpose(bi)


Hk = H_k_general(k_input)

mth_derivative_of_u = np.dot(bi_transpose, np.transpose(Hk))

l_error_est=1000

maximum_so_far = -100000

for i in range (1,l_error_est):
    p = mth_derivative_of_u.subs(t, (i - 1) / l_error_est)
    maximum_so_far = max(math.fabs(p), maximum_so_far)


M_error = maximum_so_far
print(maximum_so_far)
m_error = 1
alpha_error = 1/3
large_product = scipy.special.gamma(alpha_error)*scipy.special.gamma(m_error-alpha_error)*alpha_error*(m_error-alpha_error)*sqrt(1-pow(2,2*(alpha_error-m_error)))*pow(k_input,m_error-alpha_error)
final_error = M_error / large_product
print(final_error)





















