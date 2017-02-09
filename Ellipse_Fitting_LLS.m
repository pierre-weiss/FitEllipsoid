% function [q,A,b,c]=Ellipse_Fitting_LLS(x)
%
% Given a set of points x=(x1,..,xn), this function finds a fitting
% ellipse in 2D, by using the linear least square approach described in 
% Fitzgibbon, Andrew, Maurizio Pilu, and Robert B. Fisher. 
% "Direct least square fitting of ellipses." 
% IEEE Transactions on pattern analysis and machine intelligence 21.5 (1999): 476-480.
%
% Note : the output can be an hyperbola here, if the input points are too noisy
% or badly located on the ellipse. In that case, a warning is displayed.
%
% The output ellipse E is described implicitely by a triplet (A,b,c):
% E={x in R^2, <Ax,x> + <b,x> + c=0}
% or alternatively by vector q =(a11,a22,sqrt(2)a12,b1,b2,c).
%
% Developers : Jerome Fehrenbach & Pierre Weiss 2017.
function [q,A,b,c]=Ellipse_Fitting_LLS(x)

%% An ellipse is parameterized as a11 x^2 + a22 y^2  + 2*a12 xy +  b1 x + b2 y + c = 0
%% Vector q =(a11,a22,sqrt(2)a12,b1,b2,c)
n=size(x,2);

D=zeros(6,n);
D(1,:)=x(1,:).^2;
D(2,:)=x(2,:).^2;
D(3,:)=sqrt(2)*x(2,:).*x(1,:);
D(4,:)=x(1,:);
D(5,:)=x(2,:);
D(6,:)=1;

K=D*D';

%% The objective is now to solve min <q,Kq>, ||q||=1
[~,~,V1]=svd(K);
q=V1(:,6);
q=q/sum(q(1:2));
if q(1)<0
    q=-q;
end

if (q(2)<=0)
    disp('Beware, the LLS solution is not an ellipse');
end

A=[[q(1),q(3)/sqrt(2)];[q(3)/sqrt(2),q(2)]];
b=[q(4);q(5)];
c=q(6);
