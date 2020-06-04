import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate


R_star = 10
R_pl = 1
a = 50
i = 1.55	#np.pi/2
sample = 1000
alb = 0.52	#atmos[heric albedo


def area(theta, R):
	return (R**2) * (1/2) * (theta - np.sin(theta))


def luminosity_reflected(a, Rp, i, theta, albedo = alb):
	#phase function
	phi = (1 - np.sin(i)*np.sin(theta))/2
	return albedo * (1/4) * ((Rp / a)**2) * phi


def luminosity_thermal(a, R_pl, albedo = alb):
	return (1 - albedo) * ((R_pl/a)**2) * (1/4)


def transit_circular(a, i, theta, Rs, Rp, mode = 0, reflection = False, thermal = True):
	major_axis = a			#mode = 0 no limb darkening
	minor_axis = a * np.cos(i)	#mode = 1 linear limb darkening
					#mode = 2 quad limb darkening
	x = major_axis * np.cos(theta)	#mode = 3 3 parameter limb darkening
	y = minor_axis * np.sin(theta)	#mode = 4 4 parameter limb darkening
	d = (x**2 + y**2)**(0.5)

	L_stellar = integrate.quad(limb_darkening_integrand, 0, Rs, args = (Rs, mode))[0]

	if thermal:
		L_planet = L_stellar * luminosity_thermal(a, Rp)
	else:
		L_planet = 0

	if reflection:
		L_reflected = L_stellar * luminosity_reflected(a, Rp, i, theta)
	else:
		L_reflected = 0


	if d >= (Rs + Rp):
		A = 0
	elif d <= (Rs - Rp):
		A = np.pi * (Rp)**(2)
	else:
		theta_p = 2 * np.arccos((d**2 - (Rs**2 - Rp**2))/(2*d*Rp))
		theta_s = 2 * np.arcsin((Rp/Rs) * np.sin(theta_p / 2))
		Ap = area(theta_p, Rp)
		As = area(theta_s, Rs)
		A = Ap + As

	if mode == 0:
		intensity = 1
	elif mode == 1:
		intensity = limb_darkening_linear(min([d, Rs]), Rs)
	elif mode == 2:
		intensity = limb_darkening_quad(min([d, Rs]), Rs)
	elif mode == 3:
		intensity = limb_darkening_nlinear3(min([d, Rs]), Rs)
	elif mode == 4:
		intensity = limb_darkening_nlinear4(min([d, Rs]), Rs)
	else:
		intensity = 0

	L_blocked = A * intensity
	L_total = L_stellar + L_reflected + L_planet

	return (L_total - L_blocked)


def occultation_circular(a, i, theta, Rs, Rp, mode = 0, albedo = alb, thermal = False):
	major_axis = a
	minor_axis = a * np.cos(i)

	x = major_axis * np.cos(theta)
	y = minor_axis * np.sin(theta)
	d = (x**2 + y**2)**(0.5)

	L_stellar = integrate.quad(limb_darkening_integrand, 0, Rs, args = (Rs, mode))[0]

	if thermal:
		L_planet = luminosity_thermal(a, Rp) * L_stellar
	else:
		L_planet = 0


	if d >= (Rs + Rp):
		A = 0
	elif d <= (Rs - Rp):
		A = np.pi * (Rp)**2
	else:
		theta_p = 2 * np.arccos((d**2 - (Rs**2 - Rp**2))/(2*d*Rp))
		theta_s = 2 * np.arcsin((Rp/Rs) * np.sin(theta_p / 2))
		Ap = area(theta_p, Rp)
		As = area(theta_s, Rs)
		A = As + Ap

	L_blocked = (albedo/(4*np.pi)) * (A/(a**2)) * L_stellar
	L_reflected = L_stellar * luminosity_reflected(a, Rp, i, theta)
	L_planet = L_planet * (1 - (A/(np.pi * Rp**2)))

	return L_stellar + L_planet + max([L_reflected - L_blocked,0])



def limb_darkening_linear(d, R_s, u = 0.5):
	theta = np.arcsin(d/R_s)
	mu = np.cos(theta)
	return 1 - u * (1 - mu)


def limb_darkening_quad(d, R_s, u1 = 0.25, u2 = 0.4):
	theta = np.arcsin(d/R_s)
	mu = np.cos(theta)
	return 1 - u1 * (1 - mu) - u2 * (1 - mu)**2


def limb_darkening_nlinear3(d, R_s, u1 = 0.25, u2 = 0.15, u3 = 0.10):
	theta = np.arcsin(d/R_s)
	mu = np.cos(theta)
	return 1 - u1*(1 - mu) - u2*(1 - mu**(3/2)) - u3*(1 - mu**(2))


def limb_darkening_nlinear4(d, R_s, u1 = 0.30, u2 = 0.10, u3 = 0.07, u4 = 0.03):
	theta = np.arcsin(d/R_s)
	mu = np.cos(theta)
	return 1 - u1*(1 - mu**(1/2)) - u2*(1 - mu) - u3*(1 - mu**(3/2)) - u4*(1 - mu**(2))


def limb_darkening_integrand(d, R_s, mode):
	val = 4 * np.pi * R_s * d * (1 / (R_s**2 - d**2)**(1/2))
	if mode == 0:
		return val
	elif mode == 1:
		return val * limb_darkening_linear(d, R_s)

	elif mode == 2:
		return val * limb_darkening_quad(d, R_s)

	elif mode == 3:
		return val * limb_darkening_nlinear3(d, R_s)

	elif mode == 4:
		return val * limb_darkening_nlinear4(d, R_s)

	else:
		return 0



def star_intensity_plot():		#Taking into account the Limb Darkening effects
	#graphical representation
	d = np.arange(0, R_star, R_star/sample)
	y1 = [limb_darkening_linear(x, R_star) for x in d]
	y2 = [limb_darkening_quad(x, R_star) for x in d]
	y3 = [limb_darkening_nlinear3(x, R_star) for x in d]
	y4 = [limb_darkening_nlinear4(x, R_star) for x in d]

	plt.plot(d, y1, label = "Linear Model")
	plt.plot(d, y2, label = "Quadratic Model")
	plt.plot(d, y3, label = "Non Linear 3 parameters Model")
	plt.plot(d, y4, label = "NON Linear 4 parameters Model")
	plt.legend()
	plt.show()


	#Pictorial Representation
	x = np.arange(-R_star, R_star, 2 * R_star/sample)
	y = np.arange(-R_star, R_star, 2 * R_star/sample)
	I1 = []
	I2 = []
	I3 = []
	I4 = []
	for a in x:
		intensity1 = []
		intensity2 = []
		intensity3 = []
		intensity4 = []
		for b in y:
			dist = (a**2 + b**2)**(0.5)
			if dist > R_star:
				intensity1.append(0)
				intensity2.append(0)
				intensity3.append(0)
				intensity4.append(0)
			else:
				intensity1.append(limb_darkening_linear(dist, R_star))
				intensity2.append(limb_darkening_quad(dist, R_star))
				intensity3.append(limb_darkening_nlinear3(dist, R_star))
				intensity4.append(limb_darkening_nlinear4(dist, R_star))

		I1.append(intensity1)
		I2.append(intensity2)
		I3.append(intensity3)
		I4.append(intensity4)

	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
	ax1.imshow(I1, cmap = 'hot', interpolation = 'nearest')
	ax2.imshow(I2, cmap = 'hot', interpolation = 'nearest')
	ax3.imshow(I3, cmap = 'hot', interpolation = 'nearest')
	ax4.imshow(I4, cmap = 'hot', interpolation = 'nearest')
	ax1.set_title("Linear Limb Darkening")
	ax2.set_title("Quadratic Limb Darkening")
	ax3.set_title("Non Linear Limb Darkening 3 Parameters")
	ax4.set_title("Non Linear Limb Darkening 4 Parameters")
	plt.show()



def transit_plot():
	val = (((R_star + R_pl)/a)**2 - (np.cos(i))**2)**(0.5) * (1/np.sin(i))
	theta_1 = np.arccos(val)
	theta_2 = np.arccos(-1 * val)
	angle = np.arange(theta_1, theta_2, (theta_2 - theta_1)/sample)

	t_data_0 = [transit_circular(a, i, theta, R_star, R_pl, mode = 0, reflection = True) for theta in angle]
	t_data_1 = [transit_circular(a, i, theta, R_star, R_pl, mode = 1, reflection = True) for theta in angle]
	t_data_2 = [transit_circular(a, i, theta, R_star, R_pl, mode = 2, reflection = True) for theta in angle]
	t_data_3 = [transit_circular(a, i, theta, R_star, R_pl, mode = 3, reflection = True) for theta in angle]
	t_data_4 = [transit_circular(a, i, theta, R_star, R_pl, mode = 4, reflection = True) for theta in angle]

	plt.plot(angle, normalise(t_data_0), label = "TRANSIT with NO LD")
	plt.plot(angle, normalise(t_data_1), label = "TRANSIT linear LD")
	plt.plot(angle, normalise(t_data_2), label = "TRANSIT quadratic LD")
	plt.plot(angle, normalise(t_data_3), label = "TRANSIT 3 parameter LD")
	plt.plot(angle, normalise(t_data_4), label = "TRANSIT 4 parametr LD")
	plt.legend()
	plt.show()


def occultation_plot():
	val = (((R_star + R_pl)/a)**2 - (np.cos(i))**2)**(0.5) * (1/np.sin(i))
	theta_1 = np.arccos(val) + np.pi
	theta_2 = np.arccos(-1 * val) + np.pi
	angle = np.arange(theta_1, theta_2, (theta_2 - theta_1)/sample)

	t_data = [occultation_circular(a, i, theta, R_star, R_pl) for theta in angle]
	plt.plot(angle, normalise(t_data), label = "Occultation")
	plt.legend()
	plt.show()


def orbit_plot():
	x = np.arange(0, 2 * np.pi, (2 * np.pi)/sample)
	t_data_0 = []
	t_data_1 = []
	t_data_2 = []
	t_data_3 = []
	t_data_4 = []

	for theta in x:
		if theta <= np.pi:
			t_data_0.append(transit_circular(a, i, theta, R_star, R_pl, mode = 0, reflection = True))
			t_data_1.append(transit_circular(a, i, theta, R_star, R_pl, mode = 1, reflection = True))
			t_data_2.append(transit_circular(a, i, theta, R_star, R_pl, mode = 2, reflection = True))
			t_data_3.append(transit_circular(a, i, theta, R_star, R_pl, mode = 3, reflection = True))
			t_data_4.append(transit_circular(a, i, theta, R_star, R_pl, mode = 4, reflection = True))
		else:
			t_data_0.append(occultation_circular(a, i, theta, R_star, R_pl, mode = 0))
			t_data_1.append(occultation_circular(a, i, theta, R_star, R_pl, mode = 1))
			t_data_2.append(occultation_circular(a, i, theta, R_star, R_pl, mode = 2))
			t_data_3.append(occultation_circular(a, i, theta, R_star, R_pl, mode = 3))
			t_data_4.append(occultation_circular(a, i, theta, R_star, R_pl, mode = 4))


	plt.plot(x, normalise(t_data_0), label = "No limb darkening")
	plt.plot(x, normalise(t_data_1), label = "Linear LD")
	plt.plot(x, normalise(t_data_2), label = "Quad LD")
	plt.plot(x, normalise(t_data_3), label = "3 Params LD")
	plt.plot(x, normalise(t_data_4), label = "4 Params LD")
	plt.legend()
	plt.show()


def normalise(arr):
	max_val = max(arr)
	return [x / max_val for x in arr]


#star_intensity_plot()
transit_plot()
occultation_plot()

