import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate


R_star = 10
R_pl = 1
a = 50
i = 1.55	#np.pi/2
sample = 1000


def area(theta, R):
	return (R**2) * (1/2) * (theta - np.sin(theta))



def raw_transit(d, Rs, Rp):
	if d >= (Rs + Rp):
		return 1

	elif d <= (Rs - Rp):
		return 1 - (Rp/Rs)**2

	else:
		theta_p = 2 * np.arccos((d**2 - (Rs**2 - Rp**2))/(2*d*Rp))
		theta_s = 2 * np.arcsin((Rp/Rs) * np.sin(theta_p / 2))
		Ap = area(theta_p, Rp)
		As = area(theta_s, Rs)
		return 1 - (As + Ap)/(np.pi * Rs**2)


def transit_LD(d, Rs, Rp, mode = 0):	#mode = 0 no limb darkening
					#mode = 1 linear limb darkening
					#mode = 2 quad limb darkening
					#mode = 3 non linear limb darkening
	if d >= (Rs + Rp):
		return 1

	elif d <= (Rs - Rp):
		A = np.pi * (Rp)**(2)

	else:
		theta_p = 2 * np.arccos((d**2 - (Rs**2 - Rp**2))/(2*d*Rp))
		theta_s = 2 * np.arcsin((Rp/Rs) * np.sin(theta_p / 2))
		Ap = area(theta_p, Rp)
		As = area(theta_s, Rs)
		A = As + Ap


	if mode == 0:
		return 1 - A / (np.pi * Rs**2)
	elif mode == 1:
		intensity = limb_darkening_linear(min([d, Rs]), Rs)
	elif mode == 2:
		intensity = limb_darkening_quad(min([d, Rs]), Rs)
	elif mode == 3:
		intensity = limb_darkening_nlinear(min([d, Rs]), Rs)
	else:
		intensity = 0

	L_total = integrate.quad(limb_darkening_integrand, 0, Rs, args = (Rs, mode))[0]

	return 1 - (A/L_total)*intensity



def occultation(d, Rs, Rp):
	return 1



def ellipse(a, i, R_s = R_star, R_p = R_pl, mode = 0):
	minor_axis = a * np.cos(i)
	major_axis = a
	t_data = []

	for theta in np.arange(0, 2 * np.pi, (2 * np.pi)/sample):
		x = major_axis * np.cos(theta)
		y = minor_axis * np.sin(theta)
		d = ((x**2 + y**2)**(1/2))
		if theta <= np.pi:
			t_data.append(transit_LD(d, R_s, R_p, mode))
		else:
			t_data.append(occultation(d, R_s, R_p))

	return t_data



def limb_darkening_linear(d, R_s, u = 0.5):
	theta = np.arcsin(d/R_s)
	mu = np.cos(theta)
	return 1 - u * (1 - mu)


def limb_darkening_quad(d, R_s, u1 = 0.25, u2 = 0.25):
	theta = np.arcsin(d/R_s)
	mu = np.cos(theta)
	return 1 - u1 * (1 - mu) - u2 * (1 - mu)**2


def limb_darkening_nlinear(d, R_s, u1 = 0.5, u2 = 0.25, u3 = 0.25):
	theta = np.arcsin(d/R_s)
	mu = np.cos(theta)
	return 1 - u1*(1 - mu) - u2*(1 - mu**(3/2)) - u3*(1 - mu**(2))


def limb_darkening_integrand(d, R_s, mode):
	if mode == 1:
		return 2 * np.pi * d * limb_darkening_linear(d, R_s)

	if mode == 2:
		return 2 * np.pi * d * limb_darkening_quad(d, R_s)

	if mode == 3:
		return 2 * np.pi * d * limb_darkening_nlinear(d, R_s)
	return 0



def star_intensity_plot():		#Taking into account the Limb Darkening effects
	#graphical representation
	d = np.arange(0, R_star, R_star/sample)
	y1 = [limb_darkening_linear(x, R_star) for x in d]
	y2 = [limb_darkening_quad(x, R_star) for x in d]
	y3 = [limb_darkening_nlinear(x, R_star) for x in d]

	plt.plot(d, y1, label = "Linear Model")
	plt.plot(d, y2, label = "Quadratic Model")
	plt.plot(d, y3, label = "Non Linear Model")
	plt.legend()
	plt.show()


	#Pictorial Representation
	x = np.arange(-R_star, R_star, 2 * R_star/sample)
	y = np.arange(-R_star, R_star, 2 * R_star/sample)
	I1 = []
	I2 = []
	I3 = []
	for a in x:
		intensity1 = []
		intensity2 = []
		intensity3 = []
		for b in y:
			dist = (a**2 + b**2)**(0.5)
			if dist > R_star:
				intensity1.append(0)
				intensity2.append(0)
				intensity3.append(0)
			else:
				intensity1.append(limb_darkening_linear(dist, R_star))
				intensity2.append(limb_darkening_quad(dist, R_star))
				intensity3.append(limb_darkening_nlinear(dist, R_star))

		I1.append(intensity1)
		I2.append(intensity2)
		I3.append(intensity3)

	fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
	ax1.imshow(I1, cmap = 'hot', interpolation = 'nearest')
	ax2.imshow(I2, cmap = 'hot', interpolation = 'nearest')
	ax3.imshow(I3, cmap = 'hot', interpolation = 'nearest')
	ax1.set_title("Linear Limb Darkening")
	ax2.set_title("Quadratic Limb Darkening")
	ax3.set_title("Non Linear Limb Darkening")
	plt.show()



def transit_plot_partial():
	d = np.concatenate((np.arange(R_star + R_pl, a*np.cos(i), -(R_star + R_pl - a*np.cos(i))/sample),np.arange(a*np.cos(i), R_star + R_pl, (R_star + R_pl - a*np.cos(i))/sample)),axis = 0)
	t_data_r = [transit_LD(x, R_star, R_pl, mode = 0) for x in d]
	t_data_1 = [transit_LD(x, R_star, R_pl, mode = 1) for x in d]
	t_data_2 = [transit_LD(x, R_star, R_pl, mode = 2) for x in d]
	counter = range(len(d))
	plt.plot(counter, t_data_r, label = "RAW TRANSIT")
	plt.plot(counter, t_data_1, label = "TRANSIT LD1")
	plt.plot(counter, t_data_2, label = "TRANSIT LD2")
	plt.legend()
	plt.show()


def transit_plot_full():
	x = np.arange(0, 2 * np.pi, (2 * np.pi)/sample)
	t_data = ellipse(a, i, mode = 0)

	plt.plot(x, t_data, label = str(i) + "rad")
	plt.legend()
	plt.show()


star_intensity_plot()
transit_plot_partial()
transit_plot_full()


