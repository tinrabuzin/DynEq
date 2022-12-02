import scipy as sp
import numpy as np


def sglr_pert(sys, r):
    """Model reduction by singular perturbation.
    
    Parameters
    ----------
    sys : signal.StateSpace
        State-space model to be reduced.
    r : int
        Order of the reduced model
    
    Returns
    -------
    signal.StateSpace
        Reduced state-space model.
    """
    A11 = sys.A[:r,:r]
    A12 = sys.A[:r,r:]
    A21 = sys.A[r:,:r]
    A22 = sys.A[r:,r:]
    A22i = np.linalg.inv(A22)
    B1 = sys.B[:r]
    B2 = sys.B[r:]
    C1 = sys.C[:,:r] 
    C2 = sys.C[:,r:]
    
    Ar = A11 - A12@A22i@A21
    Br = B1 - A12@A22i@B2
    Cr = C1 - C2@A22i@A21
    Dr = sys.D - C2@A22i@B2
    return sp.signal.StateSpace(Ar, Br, Cr, Dr)

def trunc(sys, r):
    """Model reduction by truncation.
    
    Parameters
    ----------
    sys : signal.StateSpace
        State-space model to be reduced.
    r : int
        Order of the reduced model
    
    Returns
    -------
    signal.StateSpace
        Reduced state-space model.
    """
    Ar = sys.A[:r,:r]
    Br = sys.B[:r]
    Cr = sys.C[:,:r]
    Dr = sys.D
    return sp.signal.StateSpace(Ar, Br, Cr, Dr)

def grammians(sys):
    """Computes controllability and observability Grammians for the
    system.
    
    Parameters
    ----------
    sys : signal.StateSpace
        A State-space model.
    
    Returns
    -------
    ndarray
        Controllability Gramian
    ndarray
        Observarbility Gramian
    """
    P = sp.linalg.solve_lyapunov(sys.A, -sys.B@sys.B.T)
    Q = sp.linalg.solve_lyapunov(sys.A.T, -sys.C.T@sys.C)
    return P, Q

def hsv(sys):
    """Computes Hankel singular values of the system.
    
    Parameters
    ----------
    sys : signal.StateSpace
        A State-space model.
        
    Returns
    -------
    ndarray
        A vector of Hankel singular values.
    """
    P, Q = grammians(sys)
    return np.sqrt(np.linalg.eig(P@Q)[0])

def balance(sys):
    """Balances a state-space system with decreasingly ordered Hankel 
    singular values.
    
    Parameters
    ----------
    sys : signal.StateSpace
        A State-space model.
    
    Returns
    signal.StateSpace
        A balanced state-space model.
    ndarray
        Ordered Hankel singular values.
    """
    P, Q = grammians(sys)
    R = sp.linalg.cholesky(P, lower=True)
    U, s, _ = sp.linalg.svd(R.T@Q@R)
    T = R@U@np.diag(s**0.25)
    Tinv = np.linalg.inv(T)
    Ab = Tinv@sys.A@T
    Bb = Tinv@sys.B
    Cb = sys.C@T
    return sp.signal.StateSpace(Ab,Bb,Cb, sys.D), s

def bal_trunc(sys, r):
    """Model reduction by balanced truncation.
    
    Parameters
    ----------
    sys : signal.StateSpace
        State-space model to be reduced.
    r : int
        Order of the reduced model
        
    Returns
    -------
    signal.StateSpace
        Reduced state-space model.
    """
    bal_sys, s = balance(sys)
    return trunc(bal_sys, r), s

def bal_sglr_pert(sys, r):
    """Model reduction by balanced singular perturbation.
    
    Parameters
    ----------
    sys : signal.StateSpace
        State-space model to be reduced.
    r : int
        Order of the reduced model
        
    Returns
    -------
    signal.StateSpace
        Reduced state-space model.
    """
    bal_sys, s = balance(sys)
    return sglr_pert(bal_sys, r), s

def optim_red(sys, r):
    """Optimal Model Order Reduction in Hankel norm.
    
    Notes
    -----
        Not implemented because of the missing function to separate
        stable and unstable parts of the system!
    """
    raise NotImplementedError
    n = sys.A.shape[0]
    sysb, s = balance(sys)
    
    perm = np.vstack((np.hstack((np.eye(r), np.zeros((r,n-r)))),
                      np.hstack((np.zeros((n-r-1,r+1)), np.eye(n-r-1,n-r-1))),
                      np.zeros((1,n))
                    ))
    perm[n-1,r] = 1

    T = perm.T
    Aper = np.linalg.inv(T)@sysb.A@T
    Bper = np.linalg.inv(T)@sysb.B
    Cper = sysb.C@T
    hsv = np.sqrt(s)
    P1 = np.diag(hsv[hsv!=hsv[r+1]])
    Q1 = P1
    A11 = Aper[:n-1,:n-1]
    B1 = Bper[:n-1]
    C1 = Cper[:,:n-1]
    U = 1
    E1 = Q1@P1 - hsv[r]**2*np.eye(P1.shape[0])
    E1inv = np.diag(1/np.diag(E1))
    Ahat = E1inv@(hsv[r]**2*A11.T + Q1@A11@P1 - hsv[r]*C1.T@B1.T)
    Bhat = E1inv@(Q1@B1+hsv[r]*C1.T)
    Chat = C1@P1+hsv[r]*B1.T
    Dhat = sysb.D - hsv[r]
    Qstar = sp.signal.StateSpace(Ahat,Bhat,Chat,Dhat)

