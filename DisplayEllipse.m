% function DisplayEllipse(Rx,Ry,q,col)
%
% Given an ellipse given as an implicit equation 
% q(1) x^2 + q(2) y^2  + sqrt(2)*q(3) xy +  q(4) x + q(5) y + q(6) = 0,
% this function displays an ellipse in the current figure.
% 
% INPUT:
% - Rx, Ry : range in x and y coordinates for display.
% - q : 6x1 vector describing the ellipse.
% - col : 3x1 vector describing the color for display.
% OUTPUT: 
% - displays an ellipse in the current figure.
% 
% Developer: Pierre Weiss 2017.
function DisplayEllipse(Rx,Ry,q,col)

[X,Y]=meshgrid(linspace(Rx(1),Rx(2),100),linspace(Ry(1),Ry(2),100));
Z=q(1)*X.^2 + q(2)*Y.^2 + sqrt(2)*q(3)*X.*Y + q(4)*X+q(5)*Y+q(6);
contour(X,Y,Z,[0 0],'linewidth',2,'Color',col);
