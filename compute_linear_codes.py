"""
@author: Mikkel Hviid Thorn

Module to calculate a few facts about a code given its generator matrix. Many
functions are for making the calculations.
"""

import numpy as np
import matplotlib.pyplot as plt
import itertools as it
import time





"""
Functions Relating to Polynomials

The functions for calculus use representations in tuples for polynomials. 
poly_add and poly_mul are polynomial addition and multiplication were the 
result is not necessarily a minimum tuple representation.

The poly_div function calculates the remainder for a function over a finite
field with a prime number of elements, where the divider is a monic polynomial.
The result is a tuple with length equal to -1 of the divider.

Lastly the function translator translates between to representations of 
polynomials, one being the tuple and the other being the in string with a as 
variable name.
"""

def poly_add(f,g,p=2):
    """
    Calculates the sum over to polynomials in the field with p, a prime, 
    elements.

    Parameters
    ----------
    f : tuple
        tuple representation of polynomial.
    g : tuple
        tuple representation of polynomial.
    p : int
        prime number.

    Returns
    -------
    tuple
        tuple representation of polynomial.

    Example
    >>> poly_sum((1,1,0,1,1),(1,0,1,0,1,1))
    (1, 1, 0, 0, 0, 0)
    """
    deg = max(len(f),len(g))
    f, g = [0 for i in range(deg-len(f))] + list(f), [0 for i in range(deg-len(g))] + list(g)
    Sum = [(f[i] + g[i]) % p for i in range(deg)] # calculating sum
    return tuple(Sum)
    
def poly_mul(f,g,p=2):
    """
    Calculates the product over to polynomials in the field with p, a prime, 
    elements.

    Parameters
    ----------
    f : tuple
        tuple representation of polynomial.
    g : tuple
        tuple representation of polynomial.
    p : int
        prime number.

    Returns
    -------
    tuple
        tuple representation of polynomial.

    Example
    >>> poly_mul((1,1,0,1,1),(1,0,1,0,1,1))
    (0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1)
    """
    deg = max(len(f),len(g))
    f, g = [0 for i in range(deg-len(f))] + list(f), [0 for i in range(deg-len(g))] + list(g)
    Prod = [0 for i in range(2*deg)]
    for i in range(deg): # calculating product
        for j in range(deg):
            Prod[-i-j-1] = (Prod[-i-j-1] + f[-i-1]*g[-j-1]) % p
    return tuple(Prod)

def poly_div(f,D,p=2):
    """
    Calculates the remainder for one polynomial divided by a monic polynomial
    in the finite field with p, a prime, elements.

    Parameters
    ----------
    f : tuple
        tuple representation of polynomial.
    D : tuple
        tuple representation of the divider polynomial.
    p : int
        prime number.

    Returns
    -------
    tuple
        tuple representation of the remainder polynomial.
    
    Example
    >>> poly_div((1,0,1,0,1,1),(1,1,0,1,1),p=2)
    (0, 1, 1, 0)
    
    Example
    >>> poly_div((1,0,1,0,1,1),(1,1,0,1,1),p=3)
    (2, 2, 1, 2)
    """
    r = f
    while len(r) >= len(D): # keeps substracting the highest term
        temp = tuple([(-r[0]*d) % p for d in D]+[0 for i in range(len(r)-len(D))])
        r = poly_add(r,temp,p)
        n = 10*len(r)
        for i in range(len(r)): # finding the leading term
            if r[i] != 0:
                n = i
                break
        r = tuple(r[n:]) # reducing to optimal representation
    return tuple([0 for i in range(len(D)-len(r)-1)]+list(r))

def translator(f,mode,p=2,D='a^2+a+1'):
    """
    Translates between to different representations for polynomials. One is in
    text form, the other a tuple.

    Parameters
    ----------
    f : tuple or string
        tuple or string representation of polynomial.
    mode : 0 or 1
        tells which translation is used. 0 is from tuple to string, where 1 is
        from string to tuple.
    p : int
        prime number.
    D : tuple
        tuple representation of the divider polynomial.

    Returns
    -------
    tuple or string
        depending on the mode, the function returns the other representation.

    Example
    >>> translator((1,1,0,1,1),0)
    'a^4+a^3+a+1'
    
    >>> translator('a^4+3a^3+a+1',1,D='a^5')
    (1, 1, 0, 1, 1)
    
    >>> translator((4,0,65,0,1,155),0,p=3)
    '1a^5+2a^3+a+2'
    
    >>> translator('4a^5+65a^3+a+155',1,p=3,D='a^7')
    (0, 1, 0, 2, 0, 1, 2)
    """
    if type(D) == str:
        term_list = D.split('+')
        deg = int(term_list[0][-1])
    if type(D) == tuple:
        deg = len(D)-1
    
    if mode == 0:
        new_f = f'{f[-1] % p}'
        for i in range(2,len(f)+1):
            if f[len(f)-i] != 0:
                if f[len(f)-i] == 1:
                    if i == 2:
                        new_f = f'a+' + new_f
                    else:
                        new_f = f'a^{i-1}+' + new_f
                else:
                    if i == 2:
                        new_f = f'a+' + new_f
                    else:
                        new_f = f'{f[len(f)-i] % p}a^{i-1}+' + new_f
        return new_f
    
    if mode == 1:
        term_list = f.split('+')
        now_f = [0 for i in range(deg)]
        for term in term_list:
            if term[-1].isdigit() == True:
                if term.isdigit() == True:
                    now_f[deg-1] = int(term) % p
                elif term[0].isdigit() == True:
                    for i in range(len(term)):
                        if term[:i].isdigit() == True:
                            now_f[deg-int(term[-1])-1] = int(term[:i]) % p
                else:
                    now_f[deg-int(term[-1])-1] = 1
            else:
                if len(term) == 1:
                    now_f[deg-2] = 1
                else:
                    for i in range(len(term)):
                        if term[:i].isdigit() == True:
                            now_f[deg-2] = int(term[:i]) % p
        return tuple(now_f)
    




"""
Addition and Multiplication Maps for Finite Fields

First the addition and multiplication maps are given for finite fields with 
p, a prime, elements. Then two specific maps for the field with 4 elements is
given with string polynomial representation.

Lastly general summation and multiplication maps are given. Works best if the
ammount of elements are p^n with n>1, otherwise use the standard ones.
"""

def A(x,y,p=2,D=0):
    """
    Addition modulo p, a prime.

    Parameters
    ----------
    x : int
        integer.
    y : int
        integer.
    p : int
        prime number.
    D : tuple or int
        the monic irreducible polynomial in tuple form or 0. has no effect

    Returns
    -------
    int
        the sum modulo p.
    """
    return (x+y) % p

def M(x,y,p=2,D=0):
    """
    Multiplication modulo p, a prime.

    Parameters
    ----------
    x : int
        integer.
    y : int
        integer.
    p : int
        prime number.
    D : tuple or int
        the monic irreducible polynomial in tuple form or 0. has no effect

    Returns
    -------
    int
        the product modulo p.
    """
    return (x*y) % p

def A4(x,y,p=2,D=0):
    """
    Addition in the field with 4 elements with string representation 
    (modulo a^2+a+1).

    Parameters
    ----------
    x : string
        number in the field.
    y : string
        number in the field.
    p : int
        prime number. has no effect.
    D : tuple or int
        the monic irreducible polynomial in tuple form or 0. has no effect

    Returns
    -------
    string
        the sum in the field with 4 elements.
    """
    if x == '0':
        return y
    elif y == '0':
        return x
    elif x == y:
        return '0'
    elif x == '1':
        if y == 'a':
            return 'a+1'
        if y == 'a+1':
            return 'a'
    elif x == 'a':
        if y == '1':
            return 'a+1'
        elif y == 'a+1':
            return '1'
    elif x == 'a+1':
        if y == '1':
            return 'a'
        elif y == 'a':
            return '1'

def M4(x,y,p=2,D=0):
    """
    Multiplication in the field with 4 elements with string representation 
    (modulo a^2+a+1).

    Parameters
    ----------
    x : string
        number in the field.
    y : string
        number in the field.
    p : int
        prime number. has no effect.
    D : tuple or int
        the monic irreducible polynomial in tuple form or 0. has no effect

    Returns
    -------
    string
        the product in the field with 4 elements.
    """
    if x == '0' or y == '0':
        return '0'
    elif x == '1':
        return y
    elif y == '1':
        return x
    elif x == 'a' and y == 'a+1':
        return '1'
    elif x == 'a' and y == 'a':
        return 'a+1'
    elif x == 'a+1' and y == 'a':
        return '1'
    elif x == 'a+1' and y == 'a+1':
        return 'a'

def AA(x,y,p=2,D=(1,1,1)):
    """
    Addition map for the finite field that is the quotient space of the 
    polynomials over the finite field with p, a prime, elements and ideal 
    generated by the monic irreducible polynomial D.

    Parameters
    ----------
    x : tuple or int
        int if no polynomials are used. tuple representation for polynomials.
    y : tuple or int
        int if no polynomials are used. tuple representation for polynomials.
    p : int
        prime number.
    D : tuple
        the monic irreducible polynomial in tuple form.

    Returns
    -------
    tuple
        element in the quotient space.

    """
    return poly_div(poly_add(x,y,p),D,p)
    
def MM(x,y,p=2,D=(1,1,1)):
    """
    Multiplication map for the finite field that is the quotient space of the 
    polynomials over the finite field with p, a prime, elements and ideal 
    generated by the monic irreducible polynomial D.

    Parameters
    ----------
    x : tuple or int
        int if no polynomials are used. tuple representation for polynomials.
    y : tuple or int
        int if no polynomials are used. tuple representation for polynomials.
    p : int
        prime number.
    D : tuple
        the monic irreducible polynomial in tuple form.

    Returns
    -------
    tuple
        element in the quotient space.
    """
    return poly_div(poly_mul(x,y,p),D,p)



"""
Sum and Count in Lists

Two functions, one for summing elements in a list and the other for counting
the ammount of non-zero elements. The is_zero function is to help decide
boolian questions about zero.
"""

def is_zero(x):
    """
    Checks if some element is zero.

    Parameters
    ----------
    x : ???
        some element.

    Returns
    -------
    bool
        true if x can be characterized as zero, false otherwise.
    """
    if type(x) == int:
        if x == 0:
            return True
        else:
            return False
    elif type(x) == str:
        if x == '0':
            return True
        else:
            return False
    elif type(x) == tuple:
        if x == tuple([0 for i in range(len(x))]):
            return True
        else:
            return False
    elif type(x) == list:
        if x == [0 for i in range(len(x))]:
            return True
        else:
            return False
    else:
        return False

def sum_list(L,p=2,D=0,a=A):
    """
    Sums the elements in a list according to the summation s.

    Parameters
    ----------
    L : list
        list of some elements in a field.
    p : int
        prime number.
    D : tuple or int
        the monic irreducible polynomial in tuple form or 0.
    a : function
        addition map.

    Returns
    -------
    Sum : ???
        the sum of the elements in the list. type is determined by the 
        addition map.
    """
    Sum = L[0]
    for i in range(1,len(L)):
        Sum = a(Sum,L[i],p,D)
    return Sum

def count_nonzero(L):
    """
    Counts the number of non-zero elements in some list. Be warned that this is 
    only according to the representations of zero checked for in the condition.

    Parameters
    ----------
    L : list
        list of some elements in a field.

    Returns
    -------
    N : int
        number of non-zero elements in the list.
    """
    N = 0
    for l in L:
        if is_zero(l) == False:
            N = N+1
    return N





"""
Matrix Operations and Multiplication

All the functions uses list as representations of vectors and matrices. The 
function tp transposes a matrix. The functions mul_col and mul_row calculates
the matrix vector product for column and row vectors respectfully.
"""

def tp(G):
    """
    Transposes the matrix G

    Parameters
    ----------
    G : list
        a double list representation of a matrix.

    Returns
    -------
    list
        a double list representation of the transposed matrix of G.
    """
    return [[G[i][j] for i in range(len(G))] for j in range(len(G[0]))]

def mul_col(G,v,p=2,D=0,a=A,m=M):
    """
    Calculates the matrix vector product, with the vector being a column 
    vector. That is for a matrix G and vector v, the result if Gv.

    Parameters
    ----------
    G : list
        a double list representation of a matrix.
    v : list
        a list representation of a vector.
    p : int
        prime number.
    D : tuple or int
        the monic irreducible polynomial in tuple form or 0.
    a : function
        addition map.
    m : function
        multiplication map.

    Returns
    -------
    list
        a list representation of a vector.
    """
    return [sum_list([m(G[j][i],v[i],p,D) for i in range(len(v))],p,D,a) for j in range(len(G))]

def mul_row(G,v,p=2,D=0,a=A,m=M):
    """
    Calculates the matrix vector product, with the vector being a row vector. 
    That is for a matrix G and vector v, the result if vG.

    Parameters
    ----------
    G : list
        a double list representation of a matrix.
    v : list
        a list representation of a vector.
    p : int
        prime number.
    D : tuple or int
        the monic irreducible polynomial in tuple form or 0.
    a : function
        addition map.
    m : function
        multiplication map.

    Returns
    -------
    list
        a list representation of a vector.
    """
    return [sum_list([m(G[i][j],v[i],p,D) for i in range(len(v))],p,D,a) for j in range(len(G[0]))]





"""
Functions to Generate Codes

The first three functions generate spaces, with VS generating the vector space
over some finite field, subS generating some subspace characterized by a
generator matrix and dualS generating the dual code for some code 
characterized by a generator matrix.

The function min_d calculates the minimum distance for some code characterized
by a generator matrix, by calculating the minimum weight.
"""

def VS(n,p=2,D=0,S=[]):
    """
    Creates the n-dimensional space (set of tuples) over the set S. Given a 
    prime p it will automaticly generate the n-dimensional vector space over 
    the field with p elements. Given D non equal to 0 it will automaticly 
    generate the n-dimensional vector space over the field with p^deg(D) 
    elements.

    Parameters
    ----------
    n : int
        dimensions of the space.
    p : int
        prime number.
    D : tuple or int
        the monic irreducible polynomial in tuple form or 0.
    S : list
        a set of elements.

    Returns
    -------
    list
        the n-dimensional space (set of tuples) over the set S.
    """
    if S == [] and D != 0: # if D and p are given and S is not
        temp_S = it.product([i for i in range(p)], repeat=len(D)-1)
        S = [tuple(s) for s in temp_S]
    elif S == []: # if only p is given
        S = [i for i in range(p)]
    V = it.product(S, repeat=n)
    return [list(v) for v in V]

def subS(G,p=2,D=0,S=[],a=A,m=M):
    """
    Creates the code generated by the generator matrix G.

    Parameters
    ----------
    G : list
        a double list representation of a generator matrix.
    p : int
        prime number.
    D : tuple or int
        the monic irreducible polynomial in tuple form or 0.
    S : list
        a set of elements.
    a : function
        addition map.
    m : function
        multiplication map.

    Returns
    -------
    list
        the code.
    """
    V = VS(len(G),p,D,S)
    return [mul_row(G,v,p,D,a,m) for v in V]

def dualS(G,p=2,D=0,S=[],a=A,m=M):
    """
    Creates the dual code for the code generated by the generator matrix G.

    Parameters
    ----------
    G : list
        a double list representation of a generator matrix.
    p : int
        prime number.
    D : tuple or int
        the monic irreducible polynomial in tuple form or 0.
    S : list
        a set of elements.
    a : function
        addition map.
    m : function
        multiplication map.

    Returns
    -------
    DS : list
        the dual code.
    """
    V = VS(len(G[0]),p,D,S)
    DS = []
    for v in V: # checks which vectors are in the dual code, that is Gv^T=0
        temp = mul_col(G,v,p,D,a,m)
        if temp == [0 for i in range(len(mul_col(G,v,p,D,a,m)))]:
            DS.append(v)
        elif temp == ['0' for i in range(len(mul_col(G,v,p,D,a,m)))]:
            DS.append(v)
        elif type(v[0]) == tuple:
            if temp == [tuple([0 for i in range(len(v[0]))]) for i in range(len(mul_col(G,v,p,D,a,m)))]:
                DS.append(v)
    return DS

def min_d(G,p=2,D=0,S=[],a=A,m=M):
    """
    Calculates the minimum distance of the code generated by the generator 
    matrix G

    Parameters
    ----------
    G : list
        a double list representation of a generator matrix.
    p : int
        prime number.
    D : tuple or int
        the monic irreducible polynomial in tuple form or 0.
    S : list
        a set of elements.
    a : function
        addition map.
    m : function
        multiplication map.

    Returns
    -------
    d : int
        minimum distance for the code.
        Creates the parity check matrix for a Hamming code with parameter
    r.
    """
    C = subS(G,p,D,S,a,m)
    d = 10*len(G[0]) # automaticly sets a high bar for the minimum distance
    for c in C: # checks if a non-zero codeword has weight lower than the currect minimum
        temp = count_nonzero(c)
        if temp < d:
            temp_bool = [is_zero(x) for x in c]
            if all(temp_bool) != True:
                d = temp
    return d
 
def hamming(r,p=2,D=0,S=[],a=A,m=M):
    """
    Creates the parity check matrix for a Hamming code with parameter
    r.

    Parameters
    ----------
    r : positive integer
        parameter for the dimension of the columns.
    p : int
        prime number.
    D : tuple or int
        the monic irreducible polynomial in tuple form or 0.
    S : list
        a set of elements.
    a : function
        addition map.
    m : function
        multiplication map.

    Returns
    -------
    list
        a double list representation of a parity check matrix.
    """
    if S == [] and D != 0: # if D and p are given and S is not
        temp_S = it.product([i for i in range(p)], repeat=len(D)-1)
        S = [tuple(s) for s in temp_S]
    elif S == []: # if only p is given
        S = [i for i in range(p)]
        
    V = VS(r,p,D,S) # vector space for the columns
    L = [] # list of one vector from all one-dimensional subspaces
    for v in V:
        temp_bool = [is_zero(x) for x in v]
        if all(temp_bool) != True: # checks if it is non-zero
            temp_bool_two = True
            for l in L: # check that it is not a scaled vector of a vector in L
                for k in S:
                    if v == [m(k,x,p,D) for x in l]:
                        temp_bool_two = False
                        break
                if temp_bool_two == False:
                    break
            if temp_bool_two == True: # adds if the vector is from a new one-dimensional subspace
                L = L + [v]
    return tp(L) # transposes to get the parity check matrix





"""
Calculating Bounds
"""

def hamming_bound(q,n,d):
    """
    Function to calculate the Hamming bound given the paramteres q, n and d.

    Parameters
    ----------
    q : int
        the number of symbols in the alphabet.
    n : int
        length of the code.
    d : int
        the minimum distance of the code.

    Returns
    -------
    int
        the Hamming bound given the parameters.
    """
    if n >= d: # makes sure the condition is met
        t = int(np.floor((d-1)/2)) # defines t
        V = np.sum([np.math.comb(n,i)*((q-1)**i) for i in range(t+1)]) # calculates and defines the volume of a sphere with radius t
        return (q**n)/V

def singleton_bound(q,n,d):
    """
    Function to calculate the Singleton bound given the paramteres q, n and d.

    Parameters
    ----------
    q : int
        the number of symbols in the alphabet.
    n : int
        length of the code.
    d : int
        the minimum distance of the code.

    Returns
    -------
    int
        the Singleton bound given the parameters.
    """
    if n >= d: # makes sure the condition is met
        return q**(n-d+1)
    
def all_bounds(q,n,d):
    """
    Function to calculate and print all bounds given the paramteres q, n and d.

    Parameters
    ----------
    q : int
        the number of symbols in the alphabet.
    n : int
        length of the code.
    d : int
        the minimum distance of the code.

    Returns
    -------
    none
    """
    if n >= d:
        values = np.array([q,n,d])
        print(f'The parameters are q={q}, n={n} and d={d}.','\n')
        print(f'Hamming Bound: {hamming_bound(q,n,d)}')
        print(f'Singleton Bound: {singleton_bound(q,n,d)}')
