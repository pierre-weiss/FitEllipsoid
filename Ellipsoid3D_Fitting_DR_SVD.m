% function [q,CF,A,b,c]=Ellipsoid3D_Fitting_DR_SVD(x,nit)
%
% Given a set of points x=(x1,..,xn), this function finds a fitting
% ellipsoid in 3D, by using the approach proposed in the companion paper. 
% This method is not affine invariant. DR stands for Douglas-Rachford
% (Lions-Mercier would be more accurate).
%
% The output ellipsoid E is described implicitely by a triplet (A,b,c):
% E={x in R^3, <Ax,x> + <b,x> + c=0}
% or alternatively by vector q =(a11,a22,a33,sqrt(2)a12,sqrt(2)a13,sqrt(2)a23,b1,b2,b3,c).
%
% INPUT:
% - x: set of coordinates of size 2xn.
% - nit: number of iterations in Douglas-Rachford algorithm.
% OUTPUT:
% - q: (a11,a22,a33,sqrt(2)a12,sqrt(2)a13,sqrt(2)a23,b1,b2,b3,c).
% - CF: Cost function wrt to iterations.
% - A,b,c : matrix, vector, scalar describing the ellipse.
%

% Developers : Jerome Fehrenbach & Pierre Weiss 2017.
function [q,CF,A,b,c]=Ellipsoid3D_Fitting_DR_SVD(x,nit)

%% An ellipsoid is parameterized 
% a1 x^2 + a2 y^2 + a3 z^2 + a4 xy + a5 xz + a6 yz + a7 x + a8 y + a9 z + a10 = 0
% Vector q =(a11,a22,a33,sqrt(2)a12,sqrt(2)a13,sqrt(2)a23,b1,b2,b3,c)

n=size(x,2);

%% First finds the SVD of x and changes coordinates.
t=mean(x,2);
xb(1,:)=x(1,:)-t(1);
xb(2,:)=x(2,:)-t(2);
xb(3,:)=x(3,:)-t(3);
[U,S]=eig(xb*xb');
sp=max(diag(S),1e-15*ones(3,1));
P=diag(sp.^(-0.5))*U';

x=P*xb;

%
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

%% The objective is now to solve min <q,Kq>, Tr(Q)=1, Q>=0
c1=mean(x(1,:));
c2=mean(x(2,:));
c3=mean(x(3,:));
r2=var(x(1,:))+var(x(2,:))+var(x(3,:));

u=[1/3;1/3;1/3;0;0;0;-2*c1/3;-2*c2/3;-2*c3/3;(c1^2+c2^2+c3^2-r2)/3];

% And now go to the Douglas-Rachford (Lions-Mercier) iterative algorithm
gamma=10;  % Parameter in ]0,+infty[
M = gamma*K+eye(size(K));
proxf1= @(q) M\q;
proxf2= @(q) Project_on_B(q);
p=u;
CF=zeros(nit+1,1);
for k=1:nit
    q=proxf2(p);
    
    CF(k)=0.5*q'*K*q;
    
    p=p+1.0*(proxf1(2*q-p)-q);
end
q=proxf2(q);
CF(nit+1)=0.5*q'*K*q;
A2=[[q(1),q(4)/sqrt(2),q(5)/sqrt(2)];[q(4)/sqrt(2),q(2),q(6)/sqrt(2)];...
    [q(5)/sqrt(2),q(6)/sqrt(2),q(3)]];
b2=[q(7);q(8);q(9)];
c2=q(10);

%% Go back to the initial basis
A=P'*A2*P;
b=-2*A*t+P'*b2;
c=(A2*P*t)'*(P*t)-b2'*P*t+c2;

q=[A(1,1);A(2,2);A(3,3);sqrt(2)*A(2,1);sqrt(2)*A(3,1);sqrt(2)*A(3,2);b(1);b(2);b(3);c];
q=q/(A(1,1)+A(2,2)+A(3,3)); % Just a normalization to stay on the simplex

end

function q=Project_on_B(q0)
Q0=[[q0(1),q0(4)/sqrt(2),q0(5)/sqrt(2)];[q0(4)/sqrt(2),q0(2),q0(6)/sqrt(2)];...
    [q0(5)/sqrt(2),q0(6)/sqrt(2),q0(3)]];
[U,S0]=eig(Q0);
s0=diag(S0);
s=projsplx(s0);
S=diag(s);
Q=U*S*U';
q=[Q(1,1);Q(2,2);Q(3,3);sqrt(2)*Q(2,1);sqrt(2)*Q(3,1);sqrt(2)*Q(3,2);q0(7:end)];
end


