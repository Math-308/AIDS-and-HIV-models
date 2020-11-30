
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

beta1= 0.0927
beta2= 0.0164
k1= 1
k2= 0.0055
k4= 0.9
k5= 0.0000765
k6=0.1
delta1= 0.95 #mortality rate of AIDs
delta2= 1
N= 226.5 #the starting population
H0= 1
A0= 0
HA0= 0
AA0=0
D0= 0
S0= N-H0-A0-D0
# A grid of time points (in days)
t = np.linspace(0, 5, 480)

# The SIR model differential equations.
def deriv(y, t, beta, gamma):
    S, H, HA, A, AA, D = y
    dSdt = -beta1 * S * (H + A) -beta1 * S * (HA + AA)  # succeptibles
    dHdt = beta1 * S * (H + A) + beta2 * S * (HA + AA) - k1 * H + k2 * (H + A) -k4 * H + k5 * (HA + AA)  # HIV infected
    dHAdt = k4 * H - k6 * HA  # HIV infected
    dAdt = k1 * H - delta1 * A - k4 * A  # AIDS infected
    dAAdt = k4 * A - delta2 * AA  # AIDS infected
    dDdt = delta1 * A + delta2 * AA  # dead
    '''dSdt = -beta * S* (H+A)
    dHdt = beta * S* (H+A)- k1*H+ k2*(H+A)
    dAdt= k1*H- delta*A
    dDdt = delta*A'''
    return dSdt, dHdt, dHAdt, dAdt, dAAdt, dDdt

# Initial conditions vector
y0 = S0, H0, HA0, A0, AA0, D0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(beta1,beta2))
S, H, HA, A, AA, D = ret.T

# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S, 'g', alpha=0.5, lw=2, label='S')
ax.plot(t, H, 'b', alpha=0.5, lw=2, label='H')
ax.plot(t, HA, 'r', alpha=0.5, lw=2, label='HA')
ax.plot(t, A, 'cyan', alpha=0.5, lw=2, label='A')
ax.plot(t, AA, 'pink', alpha=0.5, lw=2, label='AA')
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