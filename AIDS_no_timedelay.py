
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

beta= 0.0927
k1= 1
k2= 0.0055
delta= 0.95 #mortality rate of AIDs
N= 226500000 #the starting population
H0= 1
A0= 0
D0= 0
S0= N-H0-A0-D0
# A grid of time points (in days)
t = np.linspace(0, 5, 480)

# The SIR model differential equations.
def deriv(y, t, beta, gamma):
    S, H, A, D = y
    dSdt = -beta * S* (H+A)
    dHdt = beta * S* (H+A)- k1*H+ k2*(H+A)
    dAdt= k1*H- delta*A
    dDdt = delta*A
    return dSdt, dHdt, dAdt, dDdt

# Initial conditions vector
y0 = S0, H0, A0, D0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(beta,delta))
S, H, A, D = ret.T

# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S, 'g', alpha=0.5, lw=2, label='S')
ax.plot(t, H, 'b', alpha=0.5, lw=2, label='H')
ax.plot(t, A, 'cyan', alpha=0.5, lw=2, label='A')
ax.plot(t, D, 'k', alpha=0.5, lw=2, label='D')

ax.set_xlabel('t (years)')
ax.set_ylabel('population')
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()