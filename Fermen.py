from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp

@dataclass
class paran:
  '''
  '''
  S0: float =100
  T: float = 35
  ti: float =0
  tf: float = 24*60

def k(paran):
  k=464,8*np.exp(-3985/(paran.T+273))
  return k

def DER(t, S, paran):
  dSdt=-k(paran)*S
  return dSdt

def model(paran, T=35):
  cond=[paran.S0]
  tempo=[[paran.ti], [paran.tf]]
  t_eval=np.arange(paran.ti, paran.tf, 100)
  sol=solve_ivp(DER, tempo, cond, t_eval=t_eval, args=(paran,))
  return sol.t, sol.y[0]
