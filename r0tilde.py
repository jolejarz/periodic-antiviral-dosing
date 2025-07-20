import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# This is the maximum time to run the numerical integration to.
max_time = 30

# This is the size of the time step.
h = 0.0001

# This is the total number of time steps.
N = int (max_time/h)

# This is the time elapsed.
t = np.zeros (N+1)

# This is the numerically integrated number of healthy cells versus time.
x = np.zeros (N+1)

# This is the numerically integrated viral load versus time.
y = np.zeros (N+1)

# This is the time between drug doses.
T = 2
T_steps = int(T/h)

def calc (tau):

	# L is the input rate of healthy cells..
	L = 100

	# d_x is the death rate of healthy cells.
	d_x = 0.1

	# beta_0 is the infection rate in the absence of the drug.
	beta_0 = 0.01

	# d_w is the death rate of immature infected cells.
	d_w = 0.1

	# d_y is the death rate of mature infected cells.
	d_y = 1

	# This is the fraction of time that the drug is active.
	f = 0.9

	# tau is the maturation time of the virus.
#	tau = 0.4
	#tau = 1.6
	tau_steps = int(tau/h)

	# r0tilde_predicted is the predicted value of the growth rate of the virus.
	r0tilde_predicted = -0.324
	#r0tilde_predicted = 0.123

	# This is the periodic infection rate.
	beta = np.zeros (T_steps)
	for i in range (int(T_steps*f), T_steps):
	    beta[i] = beta_0

	# t[0] is the initial time.
	t[0] = 0

	# x[0] is the initial abundance of healthy cells.
	x[0] = 1000

	# y[0] is the initial viral load.
	y[0] = 0.000001

	# This vector contains previous x values.
	x_previous = np.zeros (tau_steps)
	for i in range (tau_steps):
		x_previous[i] = x[0]

	# This vector contains previous y values.
	y_previous = np.zeros (tau_steps)
	for i in range (tau_steps):
		y_previous[i] = y[0]

	# Perform the integration.
	for i in range(N):

		t[i+1] = t[i] + h
		x[i+1] = x[i] + h * ( L - beta[i%T_steps]*y[i]*x[i] - d_x*x[i] ) 
		y[i+1] = y[i] + h * ( beta[(i-tau_steps)%T_steps]*y_previous[(i+1)%tau_steps]*x_previous[(i+1)%tau_steps]*np.exp(-d_w*tau) - d_y*y[i] ) 

		if x[i+1]<0:
			x[i+1]=0

		if y[i+1]<0:
			y[i+1]=0

		x_previous[(i+1)%tau_steps] = x[i+1]
		y_previous[(i+1)%tau_steps] = y[i+1]

fig, (ax1, ax2) = plt.subplots (1, 2, figsize = (14,5))

calc (0.4)

ax1.plot (t[int(10/h):int(30/h)], np.log10(y[int(10/h):int(30/h)]))
ax1.plot (t[int(10/h+T_steps/2):int(30/h+T_steps/2):T_steps], np.log10(y[int(10/h+T_steps/2):int(30/h+T_steps/2):T_steps]), color='red')
ax1.scatter (t[int(10/h+T_steps/2):int(30/h+T_steps/2):T_steps], np.log10(y[int(10/h+T_steps/2):int(30/h+T_steps/2):T_steps]), color='red')

slope, intercept, r, p, std_err = stats.linregress (t[int(10/h+T_steps/2):int(30/h+T_steps/2):T_steps], np.log(y[int(10/h+T_steps/2):int(30/h+T_steps/2):T_steps]))

slope_rounded = round (slope, 3)

ax1.set_xlabel('Time')
ax1.set_ylabel('Viral load (log'+'$_{10}$)')
ax1.set_title(r'$\tau$' + ' = 0.4;  ' + r'$\tilde{r}_0$' + ' = ' + str(slope_rounded))

calc (1.6)

ax2.plot (t[int(10/h):int(30/h)], np.log10(y[int(10/h):int(30/h)]))
ax2.plot (t[int(10/h+T_steps/2):int(30/h+T_steps/2):T_steps], np.log10(y[int(10/h+T_steps/2):int(30/h+T_steps/2):T_steps]), color='red')
ax2.scatter (t[int(10/h+T_steps/2):int(30/h+T_steps/2):T_steps], np.log10(y[int(10/h+T_steps/2):int(30/h+T_steps/2):T_steps]), color='red')

slope, intercept, r, p, std_err = stats.linregress (t[int(10/h+T_steps/2):int(30/h+T_steps/2):T_steps], np.log(y[int(10/h+T_steps/2):int(30/h+T_steps/2):T_steps]))

slope_rounded = round (slope, 3)

ax2.set_xlabel('Time')
ax2.set_ylabel('Viral load (log'+'$_{10}$)')
ax2.set_title(r'$\tau$' + ' = 1.6;  ' + r'$\tilde{r}_0$' + ' = ' + str(slope_rounded))

plt.show()
