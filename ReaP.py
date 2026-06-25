from dataclasses import dataclass
from scipy.integrate import solve_ivp
import numpy as np

@dataclass
class paran:
  '''
  Os parametros usados para solucionar o problema do bioreator sao definidos aqui
  Vm=Velocidade maxima de consumo (mg/L)
  Km= Constante de Michaellis (mg/l)
  Ea= Energia de ativação (j/mol)
  R= Constante universal dos gas (j/mol K)
  Tr= Temperatura de referencia (K)
  '''
  Vm: float = 5
  Km: float = 20
  Ea: float = 30000
  R: float = 8,314
  Tr: float = 298,15

def Vmax(T, paran):
  Vmaximo=paran.Vm*np.exp(-(paran.Ea/paran.R)*((1/T)-(1/paran.Tr)))
  return Vmaximo

def reator(t,S, T,paran):
  dSdt=-Vmax(T, paran)*(S/(paran.Km+S))
  return dSdt

def simu(t_points, So, T_value, param_instance):
  t_span = [t_points[0], t_points[-1]] # Span of the integration
  sol = solve_ivp(reator, t_span, [So], t_eval=t_points, args=(T_value, param_instance))
  return sol.t, sol.y[0]
