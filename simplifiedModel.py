from pylab import array, linspace, subplots
from ddeint import ddeint

# We solve the following system:
# X(t) = 1 (t < 0)
# Y(t) = 2 (t < 0)
# dX/dt = X * (1 - Y(t-d)) / 2
# dY/dt = -Y * (1 - X(t-d)) / 2

beta= 0.0927
k1= 1
k2= 0.0055
delta= 0.95 #mortality rate of AIDs
tau1= 10 #time in years
tau2= 0.75 #time in years
tau3= 3 #time delay in years
N= 226.5 #the starting population
I0= 1

D0= 0
S0= N-I0-D0
def model(Y, t, tau1, tau2, tau3):
    S, I, D = Y (t)
    S1, I1, D1 = Y(t-tau1)
    S2, I2, D2 = Y(t-tau2)
    S3, I3, D3 = Y(t-tau3)
    print('------------------------')

    dSdt= -beta * S* (I) #succeptibles
    dIdt= beta * S* (I)- +k2*(I2) - delta*I2  #infected
    dDdt= delta*I2 #dead
    return array([dSdt, dIdt, dDdt])


g = lambda t: array([S0, I0, D0]) #initial values
tt = linspace(14, 40, 8000)

fig, ax = subplots(1, figsize=(4, 4))

yy = ddeint(model, g, tt, fargs=(tau1, tau2, tau3))
ax.plot(tt, yy[:, 0], color='r', label='S')
ax.plot(tt, yy[:, 1], color='g', label='I')
ax.plot(tt, yy[:, 2], color='b', label='D')
ax.figure.savefig("figure2.jpeg")
