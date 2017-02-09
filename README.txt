This toolbox contains a set of functions to fit ellipses or ellipsoids algebraically. 
The algorithms are described precisely in 
FitEllipsoid: a user-friendly supervised segmentation plugin.
Jérôme Fehrenbach, Bastien Kovac & Pierre Weiss
Preprint, February 2017.


Detailed contents: 

Fitting ellipses in 2D: 
- function [q,CF,A,b,c]=Ellipse_Fitting_DR(x,nit)
 Approach proposed in the paper. Not affine invariant.
- function [q,CF,A,b,c]=Ellipse_Fitting_DR_SVD(x,nit)
 Approach proposed in the paper. Affine invariant.
- function [q,A,b,c]=Ellipse_Fitting_LLS(x)
 Approach of Fitzgibbons etal. Not affine invariant.
- function [q,A,b,c]=Ellipse_Fitting_LLS_SVD(x)
 Approach of Fitzgibbons etal. Affine invariant. 
- function [q,A,b,c]=Ellipse_Fitting_ALS(x,sigma)
 Approach of Kukush etal. Similar to LLS, but consistent when additional noise is added.

Fitting ellipsoids in 3D: 
- function [q,CF,A,b,c]=Ellipsoid_Fitting_DR(x,nit)
 Approach proposed in the paper. Not affine invariant.
- function [q,CF,A,b,c]=Ellipsoid_Fitting_DR_SVD(x,nit)
 Approach proposed in the paper. Affine invariant.
- function [q,A,b,c]=Ellipsoid_Fitting_LLS(x)
 Approach of Fitzgibbons etal. Not affine invariant.
- function [q,A,b,c]=Ellipsoid_Fitting_LLS_SVD(x)
 Approach of Fitzgibbons etal. Affine invariant. 

Test: 
- XP_CompareMethods_paper:
 A script to reproduce the experiments that compare the different approaches in 2D.

Additional tools:
- function DisplayEllipse(Rx,Ry,q,col)
 Displays an ellipse given as an implicit equation through a vector q.
- function DisplayEllipsoid(Rx,Ry,Rz,q,col)
 Displays an ellipsoid given as an implicit equation through a vector q.