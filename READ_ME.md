Global Variables Description:

R_star: Radius of star
R_pl: Radius of Planet
a: Orbital Distance
i: angle of inclination
sample: no of samples used in plots
alb: atmospheric albedo of the planet


Function description

area(theta, R):
    Gives the area of the segment which subtend an angle of theta in a circle of radius of R


luminosity_reflected(a, Rp, i, theta, albedo):
    Gives the fraction of luminosity reflected based on atmospheric albedo and orbital phase function for circular orbit


luminosity_thermal(a, R_pl, albedo):
    Gives the inherent luminosity of the planet based on thermal equilibrium with the star


transit_circular(a, i, theta, Rs, Rp, mode, reflection, thermal):
    Calculates the transit luminosity for the star for a circular orbit. Mode decides the model used for calculation of Limb Darkening effects.
Mode = 0    No limb darkening
Mode = 1    Linear Limb darkening model
Mode = 2    Quadratic limb darkening model
Mode = 3    Non linear limb darkening with 3 parameters
Mode = 4    Non Linear limb darkening with 4 parameters
    reflection is a boolean parameter which decides if reflected light from the planet has to be considered.
Reflection = True   Planet refelcts host star's light
Reflection = False  Planet does not reflect
    thermal is a boolean parameter which decides if the planet has luminosity of its own.
Thermal = True  Planet has its own luminosity calculated from thermal equilibrium condition with star
Thermal = False Planet does not emit 


occultation_circular(a, i, theta, Rs, Rp, mode, albedo, thermal):
    Calculates the occultation luminosity for the star for a circular orbit. Reflective property of the planet is inherently assumed.  Mode decides the model used for calculation of Limb Darkening effects.
Mode = 0    No limb darkening
Mode = 1    Linear Limb darkening model
Mode = 2    Quadratic limb darkening model
Mode = 3    Non linear limb darkening with 3 parameters
Mode = 4    Non Linear limb darkening with 4 parameters
     thermal is a boolean parameter which decides if the planet has luminosity of its own.
Thermal = True  Planet has its own luminosity calculated from thermal equilibrium condition with star
Thermal = False Planet does not emit 


normalise(arr):
    Utility function to normalize an array with respect to its peak value.

star_intensity_plot():
    plots the limb darkening effects for the above three models using the given parameters.

transit_plot():
    plots transit curve over the period of transit

occultation_plot():
    plots occultation curve over the period of occultation

orbit_plot():
     plots transit curve over the entire orbital period (transit + occultation)



