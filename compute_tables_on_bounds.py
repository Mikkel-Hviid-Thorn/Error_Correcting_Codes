"""
@author: Mikkel Hviid Thorn

Computing tables for possible values on A_q(n,d) and B_q(n,d).
"""



import numpy as np

def V(q,n,r):
    """
    Volumen of a sphere that contain words of length n and with radius r, over 
    an alphabet with q elements.

    Parameters
    ----------
    q : int
        Number of elements in the alphabet.
    n : int
        Length of the words.
    r : int
        Radius of the sphere.

    Returns
    -------
    int
        The volume of the sphere.
    """
    return np.sum([np.math.comb(n,i)*((q-1)**i) for i in range(r+1)])



# Bounds

def Hamming(q,n,d):
    """
    Computes the upper bound called the Hamming bound.

    Parameters
    ----------
    q : int
        Number of elements in the alphabet.
    n : int
        Length of the words.
    d : int
        The lower bound on the minimum distance.

    Returns
    -------
    int
        The Hamming bound for the parameters.
    """
    t = int(np.floor((d-1)/2))
    return np.floor((q**n)/V(q,n,t))

def Singleton(q,n,d):
    """
    Computes the upper bound called the Singleton bound.

    Parameters
    ----------
    q : int
        Number of elements in the alphabet.
    n : int
        Length of the words.
    d : int
        The lower bound on the minimum distance.

    Returns
    -------
    int
        The Singleton bound for the parameters.
    """
    return np.floor(q**(n-d+1))

def Gilbert(q,n,d):
    """
    Computes the lower bound called the Gilbert bound.

    Parameters
    ----------
    q : int
        Number of elements in the alphabet.
    n : int
        Length of the words.
    d : int
        The lower bound on the minimum distance.

    Returns
    -------
    int
        The Gilbert bound for the parameters.
    """
    return np.ceil((q**n)/V(q,n,d-1))

def Varshamov(q,n,d):
    """
    Computes the lower bound on linear codes called the Varshamov bound.

    Parameters
    ----------
    q : int
        Number of elements in the alphabet.
    n : int
        Length of the words.
    d : int
        The lower bound on the minimum distance.

    Returns
    -------
    int
        The Varshamov bound for the parameters.
    """
    return np.ceil(q**(n-np.ceil(np.math.log(1+V(q,n-1,d-2),q))))



# Comparing bounds

def compare_upper(q,n,d):
    """
    Compares the two upper bounds, the Hamming and Singleton bound, and uses
    results on the maximum sizes of codes to improve these.

    Parameters
    ----------
    q : int
        Number of elements in the alphabet.
    n : int
        Length of the words.
    d : int
        The lower bound on the minimum distance.

    Returns
    -------
    int
        The lowest upper bound attained. If n<d, then returns 1.
    """
    if n < d:
        return 1
    else:
        B = [q**n]
        B = B + [Hamming(q,n,d),Hamming(q,n-1,d-1),q*Hamming(q,n-1,d)]
        B = B + [Singleton(q,n,d),Singleton(q,n-1,d-1),q*Singleton(q,n-1,d)]
        return int(min(B))
    
def compare_lower(q,n,d):
    """
    Compares the two lower bounds, the Gilbert and Varshamov bound, and uses
    results on the maximum sizes of codes to improve these.

    Parameters
    ----------
    q : int
        Number of elements in the alphabet.
    n : int
        Length of the words.
    d : int
        The lower bound on the minimum distance.

    Returns
    -------
    int
        The biggest lower bound attained. If n<d, then returns 1.
    """
    if n < d:
        return 1
    else:
        B = [q]
        B = B + [Gilbert(q,n,d),Gilbert(q,n+1,d+1),Gilbert(q,n+1,d)/q]
        B = B + [Varshamov(q,n,d),Varshamov(q,n+1,d+1),Varshamov(q,n+1,d)/q]
        return int(max(B))



# Constructing tables

def table_2():
    """
    Creates a table on upper and lower bounds on the maximum sizes of codes in
    the case q=2. The table runs over parameters 4<n<21, d=4,6,8,10.
    """
    
    # Storing information
    T_upper = [[0 for i in range(4)] for j in range(16)]
    T_lower = [[0 for i in range(4)] for j in range(16)]
    
    # Calculating bounds
    for d_half in range(2,6):
        for n in range(5,21):
            T_upper[n-5][d_half-2] = compare_upper(2,n,2*d_half)
            T_lower[n-5][d_half-2] = compare_lower(2,n,2*d_half)

    # Printing table
    print('n   d=4   d=6   d=8   d=10\n')
    for n in range(5,21):
        print(f'{n}   {T_lower[n-5][0]}-{T_upper[n-5][0]}   {T_lower[n-5][1]}-{T_upper[n-5][1]}   {T_lower[n-5][2]}-{T_upper[n-5][2]}   {T_lower[n-5][3]}-{T_upper[n-5][3]}')

def table_q(q):
    """
    Creates a table on upper and lower bounds on the maximum sizes of codes. 
    The table runs over parameters 4<n<16, 3<d<8.
    """
    
    # Storing information
    T_upper = [[0 for i in range(4)] for j in range(11)]
    T_lower = [[0 for i in range(4)] for j in range(11)]
    
    # Calculating bounds
    for d in range(4,8):
        for n in range(5,16):
            T_upper[n-5][d-4] = compare_upper(q,n,d)
            T_lower[n-5][d-4] = compare_lower(q,n,d)

    # Printing table
    print('n   d=4   d=5   d=6   d=7\n')
    for n in range(5,16):
        print(f'{n}   {T_lower[n-5][0]}-{T_upper[n-5][0]}   {T_lower[n-5][1]}-{T_upper[n-5][1]}   {T_lower[n-5][2]}-{T_upper[n-5][2]}   {T_lower[n-5][3]}-{T_upper[n-5][3]}')
