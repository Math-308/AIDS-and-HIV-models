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
tau1= 0#10 #time in years
tau2= 0#0.75 #time in years
tau3= 0#3 #time delay in years
N= 226500000 #the starting population
H0= 1
A0= 0
D0= 0
S0= N-H0-A0-D0
def model(Y, t, tau1, tau2, tau3):
    S, H, A, D = Y (t)
    S1, H1, A1, D1 = Y(t-tau1)
    S2, H2, A2, D2 = Y(t-tau2)
    S3, H3, A3, D3 = Y(t-tau3)
    print('------------------------')
    print(S, H, A, D)
    print(S1, H1, A1, D1)
    print(S2, H2, A2, D2)
    print(S3, H3, A3, D3)
    print('------------------------')
    dSdt= -beta * S* (H+A) #succeptibles
    dHdt= beta * S* (H+A)- k1*H1 +k2*(H2+A2) #HIV infected
    dAdt=k1*H1 - delta*A3 #AIDS infected
    dDdt= delta*A3 #dead
    return array([dSdt, dHdt, dAdt, dDdt])


g = lambda t: array([S0, H0, A0, D0]) #initial values
tt = linspace(11, 40, 96000)

fig, ax = subplots(1, figsize=(4, 4))

yy = ddeint(model, g, tt, fargs=(tau1, tau2, tau3))
ax.plot(tt, yy[:, 0], color='r', label='S')
ax.plot(tt, yy[:, 1], color='g', label='H')
ax.plot(tt, yy[:, 2], color='b', label='A')
ax.plot(tt, yy[:, 3], color='purple', label='D')
ax.figure.savefig("figure1.jpeg")