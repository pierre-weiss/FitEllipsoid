% function [q,A,b,c]=Ellipse_Fitting_LLS_SVD(x)
%
% Given a set of points x=(x1,..,xn), this function finds a fitting
% ellipse in 2D, by:
% i) changing the coordinate system to get an affine invariant result
% ii) using the linear least square approach described in 
% Fitzgibbon, Andrew, Maurizio Pilu, and Robert B. Fisher. 
% "Direct least square fitting of ellipses." 
% IEEE Transactions on pattern analysis and machine intelligence 21.5 (1999): 476-480.
%
% Note : the output can be an hyperbola here, if the input points are too noisy
% or badly located on the ellipse. In that case, a warning is displayed.
%
% The output ellipse E is described implicitly by a triplet (A,b,c):
% {x in R^2, <Ax,x> + <b,x> + c=0}
% or alternatively by vector q =(a11,a22,sqrt(2)a12,b1,b2,c).
%
% Developers : Jerome Fehrenbach & Pierre Weiss 2017.
function [q,A,b,c]=Ellipse_Fitting_LLS_SVD(x)

%% First finds the SVD of x and changes coordinates.
t=mean(x,2);
xb(1,:)=x(1,:)-t(1);
xb(2,:)=x(2,:)-t(2);
[U,S]=eig(xb*xb');
sp=max(diag(S),1e-10*ones(2,1));
P=diag(sp.^(-0.5))*U';

x=P*xb;

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
    disp('Beware, the LLS-SVD solution is not an ellipse');
end

%% Goes back to the initial basis
A2=[[q(1),q(3)/sqrt(2)];[q(3)/sqrt(2),q(2)]];
b2=[q(4);q(5)];
c2=q(6);

A=P'*A2*P;
b=-2*A*t+P'*b2;
c=(A2*P*t)'*(P*t)-b2'*P*t+c2;

q=[A(1,1);A(2,2);sqrt(2)*A(1,2);b(1);b(2);c];
q=q/(A(1,1)+A(2,2)); % Just a normalization to stay on the simplex