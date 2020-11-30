from pylab import array, linspace, subplots
from ddeint import ddeint

# We solve the following system:
# X(t) = 1 (t < 0)
# Y(t) = 2 (t < 0)
# dX/dt = X * (1 - Y(t-d)) / 2
# dY/dt = -Y * (1 - X(t-d)) / 2

beta1= 0.0927
beta2= 0.0164
k1= 1
k2= 0.0055
k4= 0.9
k5= 0.0000765
k6=0.1
delta1= 0.95 #mortality rate of AIDs
delta2= 1
tau1= 10 #time in years
tau2= 0.75 #time in years
tau3= 35 #time delay in years
tau4= 3
tau5= 10
N= 226.5 #the starting population
H0= 1
A0= 0
HA0= 0
AA0=0
D0= 0
S0= N-H0-A0-D0
def model(Y, t, tau1, tau2, tau3, tau4, tau5):
    S, H, HA, A, AA, D = Y (t)
    S1, H1, HA1, A1, AA1, D1 = Y (t-tau1)
    S2, H2, HA2, A2, AA2, D2 = Y (t-tau2)
    S3, H3, HA3, A3, AA3, D3 = Y (t-tau3)
    S4, H4, HA4, A4, AA4, D4 = Y (t-tau4)
    S5, H5, HA5, A5, AA5, D5 = Y (t-tau5)

    dSdt= -beta1 * S* (H+A)- -beta1 * S* (HA+AA) #succeptibles
    dHdt= beta1 * S* (H+A)+ beta1 * S* (HA+AA)- k1*H1 +k2*(H2+A2)- k4*H+ k5*(HA2+AA2) #HIV infected
    dHAdt = k4*H- k6*HA3 # HIV infected
    dAdt=k1*H1 - delta1*A4- k4*A #AIDS infected
    dAAdt=k4*A - delta2*AA5#AIDS infected
    dDdt= delta1*A4+ delta2*AA5 #dead
    return array([dSdt, dHdt, dHAdt, dAdt,dAAdt, dDdt])


g = lambda t: array([S0, H0, HA0, A0, AA0, D0]) #initial values
tt = linspace(14, 40, 480000)

fig, ax = subplots(1, figsize=(4, 4))

yy = ddeint(model, g, tt, fargs=(tau1, tau2, tau3))
ax.plot(tt, yy[:, 0], color='r', label='S')
ax.plot(tt, yy[:, 1], color='g', label='H')
ax.plot(tt, yy[:, 2], color='b', label='HA')
ax.plot(tt, yy[:, 3], color='purple', label='A')
ax.plot(tt, yy[:, 4], color='b', label='AA')
ax.plot(tt, yy[:, 5], color='b', label='D')


ax.figure.savefig("figure_3.jpeg")