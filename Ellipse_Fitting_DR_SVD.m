% function [q,CF,A,b,c]=Ellipse_Fitting_DR_SVD(x,nit)
%
% Given a set of points x=(x1,..,xn), this function finds a fitting
% ellipse in 2D, by using the approach proposed in the companion paper.
% This method is affine invariant: it starts by changing the coordinate
% using the SVD.  DR stands for Douglas-Rachford (Lions-Mercier would be more accurate).
%
% The output ellipse E is described implicitely by a triplet (A,b,c):
% E={x in R^2, <Ax,x> + <b,x> + c=0}
% or alternatively by vector q =(a11,a22,sqrt(2)a12,b1,b2,c).
%
% INPUT:
% - x: set of coordinates of size 2xn.
% - nit: number of iterations in Douglas-Rachford algorithm.
% OUTPUT:
% - q: (a11,a22,sqrt(2)a12,b1,b2,c).
% - CF: Cost function wrt to iterations.
% - A,b,c : matrix, vector, scalar describing the ellipse.
%
% Example:
% t=linspace(0,2*pi,100)';
% x=[];
% x(1,:)=2*cos(t)+0.2*randn(100,1);
% x(2,:)=sin(t)+0.2*randn(100,1);
% [q,CF,A,b,c]=Ellipse_Fitting_DR_SVD(x,100);
% figure(1);hold off;
% plot(x(1,:),x(2,:),'k*');
% hold on;
% DisplayEllipse([-3,3],[-3,3],q,[1 0 0]);
% axis equal; hold off;
% legend('Input points','Fitted ellipse')
%
% Developers : Jerome Fehrenbach & Pierre Weiss 2017.
function [q,CF]=Ellipse_Fitting_DR_SVD(x,nit)

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

%% The objective is now to solve min <q,Kq>, Tr(Q)=1, Q>=0
c1=0;
c2=0;
r2=1;

u=[0.5;0.5;0;-c1;-c2;(c1^2+c2^2-r2)/2];
% And now go to the Douglas-Rachford iterative algorithm
gamma=100;  % Parameter in ]0,+infty[
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

%% Goes back to the initial basis
A2=[[q(1),q(3)/sqrt(2)];[q(3)/sqrt(2),q(2)]];
b2=[q(4);q(5)];
c2=q(6);

A=P'*A2*P;
b=-2*A*t+P'*b2;
c=(A2*P*t)'*(P*t)-b2'*P*t+c2;

q=[A(1,1);A(2,2);sqrt(2)*A(1,2);b(1);b(2);c];
q=q/(A(1,1)+A(2,2)); % Just a normalization to stay on the simplex
end

function q=Project_on_B(q0)
Q0=[[q0(1),q0(3)/sqrt(2)];[q0(3)/sqrt(2),q0(2)]];
[U,S0]=eig(Q0);
s0=diag(S0);
s=projsplx(s0);
S=diag(s);
Q=U*S*U';
q=[Q(1,1);Q(2,2);sqrt(2)*Q(1,2);q0(4:end)];
end
