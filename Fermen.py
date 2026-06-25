from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp

@dataclass
class paran:
  '''
  '''
  S0: float =100
  T: float = 35

def k(paran):
  k=464,8*np.exp(-3985/(paran.T+273))
  return k

def DER(t, S, paran):
  dSdt=-k(paran)*S
  return dSdt

def model(t, paran.S0, paran, T=35):
  cond=[paran.S0]
  T= paran.T
  t_eval=np.arange(0, t, 1000)
  sol=solve_ivp(DER, [t[0],t[-1]], cond, t_eval=t_eval, args=(paran))
  return sol.t, sol.y[0]
