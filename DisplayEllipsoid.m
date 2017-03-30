% function DisplayEllipsoid(Rx,Ry,Rz,q,col)
%
% Given an ellipse given as an implicit equation
% q(1) x^2 + q(2) y^2  + q(3) z^2 + sqrt(2)*q(4) xy + sqrt(2)*q(5) xz  +  sqrt(2)*q(6) yz  +  q(7) x + q(8) y + q(9) z + q(10) = 0,
% this function displays an ellipse in the current figure.
%
% INPUT:
% - Rx, Ry, Rz : range in x,y,z coordinates for display.
% - q : 10x1 vector describing the ellipse.
% - col : 3x1 vector describing the color for display.
% OUTPUT:
% - displays an ellipsoid in the current figure.
%
% Developer: Pierre Weiss 2017.
function DisplayEllipsoid(Rx,Ry,Rz,q,col)


[X,Y,Z]=meshgrid(linspace(Rx(1),Rx(2),200),linspace(Ry(1),Ry(2),200),linspace(Rz(1),Rz(2),200));
%[X,Y,Z]=meshgrid(linspace(Rx(1),Rx(2),10),linspace(Ry(1),Ry(2),10),linspace(Rz(1),Rz(2),10));


V=q(1)*X.^2+q(2)*Y.^2+q(3)*Z.^2+sqrt(2)*q(4)*X.*Y+sqrt(2)*q(5)*X.*Z+sqrt(2)*q(6)*Y.*Z+...
    q(7)*X+q(8)*Y+q(9)*Z+q(10);
p=patch(isosurface(X,Y,Z,V,0));
isonormals(X,Y,Z,V,p);
set(p,'FaceColor',col,'EdgeColor','none');
set(p,'FaceAlpha',0.5);
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud
hold off;
