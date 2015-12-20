# RAJZEL
RAJZEL

2D parallel first-arrival traveltime modelling and inversion code using an Eikonal solver and the adjoint state method, 
which I developed together with Denise De Nil.

The Eikonal solver is based on a fast sweeping method (Zhao, 2004) and the adjoint problem solved according to Leung & Qian (2006), Taillandier et al. (2009) and Bretaudeau et al. (2014). Intial linear 1D gradient models for the First-Arrival Traveltime Tomography (FATT) can be estimate by a grid search approach. Possible optimization methods are quasi Newton l-BFGS or non-linear PCG (Nocedal & Wright, 2006) in combination with a simple line search algorithm satisfying the weak Wolfe conditions (Skajaa, 2010). 

RAJZEL is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 2.0 of the License only.

RAJZEL is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License in LICENSE.md for more details.

If you show modelling/inversion results in a paper or presentation please 
give a reference to the following papers:

--

Daniel Koehn
