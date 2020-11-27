import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
# modeling AIDS in the US in 1980 no one is treated
# Total population, N.
N = 226.5* 1000000
# Information regarding each member of the population
#status = [] * N
# Initial number of infected and recovered individuals, I0 and R0.
H0, A0, D0 = 1, 0, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - H0 - A0- D0
# infection rate, beta(in years).
beta= 0.0927
#t is time measured in years
K1= 1 # likelihood HIV will progress to AIDS
K2= 0.0055 #likelihood a baby is born that has aids
delta= 0.95 #likelihood you die from AIDS
tau1= 10 # time that you will remain in HIV stage
tau2= 0.75 #gestation period for baby
tau3= 3 #time you have AIDS before you die
# A grid of time points (in years) from 1980 to 2020

t = np.linspace(0, 40, 40)

# The SIR model differential equations.
def deriv(y, t, N, beta, gamma):
    S, H, A, D = y
    dSdt = -beta * S * (H+A)
    dHdt = beta * S * (H+A) - K1* H(t-tau1)+ K2(H+A)*(t-tau2)
    dAdt = K1* H(t-tau1)- delta*A*(t-tau3)
    dDdt= delta*A*(t-tau3)
    return dSdt, dIdt, dRdt

# Initial conditions vector
y0 = S0, H0, A0, D0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, delta))
S, H, A, D= ret.T

# Plot the data on three separate curves for S(t), I(t) and, R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S, 'g', alpha=0.5, lw=2, label='S')
ax.plot(t, H, 'b', alpha=0.5, lw=2, label='H')
ax.plot(t, A, 'r', alpha=0.5, lw=2, label='A')
ax.plot(t, D, 'bl', alpha=0.5, lw=2, label='R')

ax.set_xlabel('t (years)')
ax.set_ylabel('f(t)')
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()