Global Variables Description:

R_star: Radius of star
R_pl: Radius of Planet
a: Orbital Distance
i: angle of inclination
sample: no of samples used in plots


Function description

area(theta, R):
    Gives the area of the segment which subtend an angle of theta in a circle of radius of R


raw_transit(d, Rs, Rp):
    Calculates the overlapping area of the planet and the host star. Uniform luminosity is assumed, leading the resultant relative luminosity to be (1 - overlap area/total area)


transit_LD(d, Rs, Rp, mode):
    Calculates the relative luminosity assuming limb darkening effects. Mode defines the model used for limb darkening of star.
Mode = 0    No limb darkening
Mode = 1    Linear Limb darkening model
Mode = 2    Quadratic limb darkening model
Mode = 3    Non linear limb darkening with 3 parameters


ellipse(a, i, R_s, R_p):
    Calculates the elliptic projection of circular orbit and returns the transit curve.


star_intensity_plot():
    plots the limb darkening effects for the above three models using the given parameters.

transit_plot_partial():
    plots transit curve over the period of transit (occultation is not considered)

transit_plot_fulll():
     plots transit curve over the entire orbital period (occultation is not considered)


