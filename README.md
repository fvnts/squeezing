# squeezing
Supplementary Materials for arXiv:1704.05840



- - - - - - 

last revision: 10.April.2017
last retrieved: 7.October.2017

developed by jes'us fuentes
email: jfuentes [at] fis.cinvestav.mx

- - - - - - 


- - - - - - 
1
- - - - - -

use genBeta.nb in Mathematica so as to 
modify the parameters of a beta-function

use genBeta.m in Matlab so as to generate
the suitable parameters that allow to have a 
proper beta-function for amplification/squeezing.
this means that beta(t)>=0 for all t in [-pi/2,pi/2]
the present code is implemented for the particular 
case of beta0 = 0.

- - - - - - 
2
- - - - - -

use evolutionBetaBeta.m in Matlab to generate the 
data files that provide a smooth evolution of a 
successive application of two non-constant beta-fields.

use evolutionBetaConst.m in Matlab to generate the 
data files that provide a smooth evolution of a 
successive application of a non-constant beta-field 
and a constant beta-field. this case is obviously 
subject under the condition of beta0 > 0.

use topelitz.nb in Mathematica to obtain the maps
of the smooth evolution, either non-constants or 
non-constant + constant.

- - - - - - 
3
- - - - - -

similarly that in 2.

use shadowsBetaBeta.m or shadowsBetaConst.m in Matlab 
to generate the data files of the shadows around a given
trajectory of the field.

- - - - - - 
4
- - - - - -

use the files dropped in Floquet to
obtain the Ince–Strutt diagram related



