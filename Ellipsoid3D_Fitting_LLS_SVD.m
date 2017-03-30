% function [q,A,b,c]=Ellipsoid3D_Fitting_LLS_SVD(x)
%
% Given a set of points x=(x1,..,xn), this function finds a fitting
% ellipsoid in 3D, by:
% i) changing the coordinate system to get an affine invariant result
% ii) using the linear least square approach described in 
% Fitzgibbon, Andrew, Maurizio Pilu, and Robert B. Fisher. 
% "Direct least square fitting of ellipses." 
% IEEE Transactions on pattern analysis and machine intelligence 21.5 (1999): 476-480.
%
% Note : the output can be an hyperbola here, if the input points are too noisy
% or badly located on the ellipse. In that case, a warning is displayed.
%
% The output ellipsoid E is described implicitly by a triplet (A,b,c):
% {x in R^3, <Ax,x> + <b,x> + c=0}
% or alternatively by vector q =(a11,a22,a33,sqrt(2)a12,sqrt(2)a13,sqrt(2)a23,b1,b2,b3,c)
%
% Developers : Jerome Fehrenbach & Pierre Weiss 2017.
function [q,A,b,c]=Ellipsoid3D_Fitting_LLS_SVD(x)

%% First finds the SVD of x and changes coordinates.
t=mean(x,2);
xb(1,:)=x(1,:)-t(1);
xb(2,:)=x(2,:)-t(2);
xb(3,:)=x(3,:)-t(3);
[U,S]=eig(xb*xb');
sp=max(diag(S),1e-15*ones(3,1));
P=diag(sp.^(-0.5))*U';

x=P*xb;

%% % An ellipsoid is parameterized 
% a1 x^2 + a2 y^2 + a3 z^2 + a4 xy + a5 xz + a6 yz + a7 x + a8 y + a9 z + a10 = 0
% Vector q =(a11,a22,a33,sqrt(2)a12,sqrt(2)a13,sqrt(2)a23,b1,b2,b3,c)

n=size(x,2);

D=zeros(10,n);
D(1,:)=x(1,:).^2;
D(2,:)=x(2,:).^2;
D(3,:)=x(3,:).^2;
D(4,:)=sqrt(2)*x(1,:).*x(2,:);
D(5,:)=sqrt(2)*x(1,:).*x(3,:);
D(6,:)=sqrt(2)*x(2,:).*x(3,:);
D(7,:)=x(1,:);
D(8,:)=x(2,:);
D(9,:)=x(3,:);
D(10,:)=1;

K=D*D';

%% The objective is now to solve min <q,Kq>, ||q||=1
[~,~,V1]=svd(K);
q=V1(:,10);
q=q/sum(q(1:3));
if q(1)<0
    q=-q;
end

%% Go back to the initial basis
A2=[[q(1),q(4)/sqrt(2),q(5)/sqrt(2)];[q(4)/sqrt(2),q(2),q(6)/sqrt(2)];...
    [q(5)/sqrt(2),q(6)/sqrt(2),q(3)]];
b2=[q(7);q(8);q(9)];
c2=q(10);

A=P'*A2*P;
b=-2*A*t+P'*b2;
c=(A2*P*t)'*(P*t)-b2'*P*t+c2;

if (min(eig(A))<=0)
    disp('Beware, the LLS-SVD solution is not an ellipse');
end

q=[A(1,1);A(2,2);A(3,3);sqrt(2)*A(2,1);sqrt(2)*A(3,1);sqrt(2)*A(3,2);b(1);b(2);b(3);c];
q=q/(A(1,1)+A(2,2)+A(3,3)); % Just a normalization to stay on the simplex


