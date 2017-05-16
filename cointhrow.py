"""
a simple half-integer verlet solver for projectiles.  I've used a coin as
an example.  Includes headwind.

"""

import numpy as np
from scipy.constants import pi, g
#some initialisation of variables
rho = 1.2 #kgm^-3
C = 1.1 #assumes a thin disc perpendicular to travel direction
V_X0 = 25.0 #reasonable for a human.
V_Y0 = 0.0 #assuming a flat trajectory
Y0 = 1.8 #reasonable starting height off ground
Y_MIN = 0.2 #below this, it's missing any normal car
del_t = 0.001 #reasonable
hdwind = 1.0 #approx. 20mph in ms^-1

class Coin():
	m = None #in kgs
	A = None # in m^2
	D = None # in m
	Cd = None
	gamma = None # =F/m(v), saves space/calculation
	v_x = V_X0
	v_y = V_Y0
	v = np.sqrt((v_x**2) + (v_y**2))
	x = 0.0
	y = Y0

	def __init__(self, denomination):
		#all from royalmint.com
		dimensions = {
			"1p"   : (3.56e-03, 20.3e-03),
			"2p"   : (7.12e-03, 25.9e-03),
			"50p"  : (8.0e-03, 27.3e-03),
			"pound": (9.5e-03, 22.5e-03)
		}
		self.m, self.D = dimensions.get(denomination)
		self.A = pi * (self.D ** 2) / 4.0
		self.Cd = rho * self.A * C * 0.5
		self.gamma = self.Cd / self.m

	def update_v(self, del_t):
		a_x = - (self.gamma * (self.v * (self.v_x + hdwind)))
		a_y = -g - (self.gamma * (self.v * self.v_y))
		self.v_x += (a_x * del_t)
		self.v_y += (a_y * del_t)
		self.v = np.sqrt((self.v_x**2) + (self.v_y**2))

	def update_r(self, del_t):
		self.x += self.v_x * del_t
		self.y += self.v_y * del_t

	def energy(self):
		return 0.5 * self.m * self.v**2

coin = Coin("pound")
print("Time\tV(x) (ms-1)\tDistance (m)\tHeight (m)\tEnergy (J)")
ticktock = True
for t in np.arange(0.001, 1.0, (del_t/2.0)):
	#use half-integer (verlet-like) method
	if (ticktock):
		coin.update_v(del_t)
	else:
		coin.update_r(del_t)
	ticktock = not ticktock
	print("{}\t{}\t{}\t{}\t{}".format(t, coin.v_x, coin.x, coin.y, coin.energy()))
	if (coin.y <= Y_MIN) or (coin.v_x <= 0.0):
		break
