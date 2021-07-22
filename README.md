# tov

In this work, I'm solving Tolman-Oppenheimer-Volkoff equations to construct a static non-rotating neutron star with a description of equation of state which could be either piecewise polytropic or, provided as table. 

Imported from natj <nattila.joonas@gmail.com>; https://github.com/natj/tov on **7.Dec.2017** into http://gitlab.icts.res.in/rahul.kashyap/_kilonovae_standardization after which these additions have been made: 
1. surface finding algorithm
2. EOS table format with interpolation
3. calculation of tidal deformability 
4. test against results from https://bitbucket.org/bernuzzi/tov/src/master/ 

#rahul/21.Jan.2019: imported from kilonovastandardization project -- https://ui.adsabs.harvard.edu/abs/2019ApJ...886L..19K/abstract 

Code 
It is still very slow but the results of non-rotating NS sequences matches up to less than 10% against the result from https://bitbucket.org/bernuzzi/tov/src/master/ which I'd recommend for high precision data. 
## To Do:
* add thermodynamimc consistent interpolation of EOS tables. 
* add numerical error in maximum mass and interpolation error. 
* Include rotation profile to the neutron stars
* Make code faster
* Make julia code; 
