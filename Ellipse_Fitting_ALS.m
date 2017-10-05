% function [q,A,b,c]=Ellipse_Fitting_ALS(x,sigma)
% 
% Given a set of points x=(x1,..,xn), this function finds a fitting
% ellipse in 2D, by using the Adjusted Least Square approach described in 
%
% Markovsky, Ivan, Alexander Kukush, and Sabine Van Huffel. 
% "Consistent least squares fitting of ellipsoids." 
% Numerische Mathematik 98.1 (2004): 177-194.
%
% The extra parameter sigma is the standard deviation of the noise. This method is consistent in the sense that the expectation of the retrieved ellipse is equal to the input ellipse. 
%
% Note : the output can be an hyperbola here, if the input points are too noisy
% or badly located on the ellipse. In that case, a warning is displayed.
%
% The output ellipse E is described implicitely by a triplet (A,b,c):
% E={x in R^2, <Ax,x> + <b,x> + c=0}
% or alternatively by vector q =(a11,a22,sqrt(2)a12,b1,b2,c).
%
% Developers : Jerome Fehrenbach & Pierre Weiss 2017.
function q=Ellipse_Fitting_ALS(x,sigma)

%% An ellipse is parameterized as q(1) x^2 + sqrt(2) * q(2) xy + q(3) y^2 + q(4) x + q(5) y + q(6) = 0
n=size(x,2);

D=zeros(6,n);
D(1,:)=x(1,:).^2;
D(2,:)=2*x(2,:).*x(1,:);
D(3,:)=x(2,:).^2;
D(4,:)=x(1,:);
D(5,:)=x(2,:);
D(6,:)=1;

K=D*D';
x1=x(1,:);
x2=x(2,:);

Delta=sigma^2*...
[[3*n*sigma^2-6*sum(x1.^2),-6*sum(x1.*x2)                ,n*sigma^2-sum(x1.^2+x2.^2),-3*sum(x1),-sum(x2)  ,-n];...
 [0                       ,4*n*sigma^2-4*sum(x1.^2+x2.^2),-6*sum(x1.*x2)            ,-2*sum(x2),-2*sum(x1), 0];...
 [0                       ,0                             ,3*n*sigma^2-6*sum(x2.^2)  ,-sum(x1)  ,-3*sum(x2),-n];...
 [0                       ,0                             ,0                         ,-n        ,0         , 0];...
 [0                       ,0                             ,0                         ,0         ,-n        , 0];...
 [0                       ,0                             ,0                         ,0         ,0         , 0]];
Delta=Delta+Delta'-diag(diag(Delta));

[V,D]=svd(K+Delta);
q=V(:,end);

q(3)=sqrt(2)*q(3);
A=[[q(1),q(3)/sqrt(2)];[q(3)/sqrt(2),q(2)]];
b=[q(4);q(5)];
c=q(6);

