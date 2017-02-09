% function [q,CF,A,b,c]=Ellipse_Fitting_DR_Kukush(x,nit,sigma)
%
% Given a set of points x=(x1,..,xn), this function finds a fitting
% ellipse in 2D, by using the approach proposed in 
%
% Markovsky, Ivan, Alexander Kukush, and Sabine Van Huffel. 
% "Consistent least squares fitting of ellipsoids." 
% Numerische Mathematik 98.1 (2004): 177-194.
%
% This method is not affine invariant, but it is consistent under additional white
% Gaussian noise of variance sigma.
%
% The output ellipse E is described implicitely by a triplet (A,b,c):
% E={x in R^2, <Ax,x> + <b,x> + c=0}
% or alternatively by vector q =(a11,a22,sqrt(2)a12,b1,b2,c).
%
% Developers : Jerome Fehrenbach & Pierre Weiss 2017.
function [q,CF,A,b,c]=Ellipse_Fitting_DR_Kukush(x,nit,sigma)

%% An ellipse is parameterized as a11 x^2 + a22 y^2  + 2*a12 xy +  b1 x + b2 y + c = 0
%% Vector q =(a11,a22,a12,b1,b2,c)
n=size(x,2);

D=zeros(6,n);
D(1,:)=x(1,:).^2;
D(2,:)=x(2,:).^2;
D(3,:)=2*x(2,:).*x(1,:);
D(4,:)=x(1,:);
D(5,:)=x(2,:);
D(6,:)=1;

K=D*D';

%% Bias correction
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
K=K+Delta;

%% The objective is now to solve min <q,Kq>, Tr(Q)=1, Q>=0
c1=mean(x(1,:));
c2=mean(x(2,:));
r2=(var(x(1,:))+var(x(2,:)));

u=[0.5;0;0.5;-c1;-c2;(c1^2+c2^2-r2)/2];
% And now go to the Lions-Mercier iterative algorithm
gamma=1;  % Parameter in ]0,+infty[
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

A=[[q(1),q(3)];[q(3),q(2)]];
b=[q(4);q(5)];
c=q(6);

q(3)=sqrt(2)*q(3); % In Kukush's method's, the coefficient vector q(3)=a12.

end

function q=Project_on_B(q0)
Q0=[[q0(1),q0(3)];[q0(3),q0(2)]];
[U,S0]=eig(Q0);
s0=diag(S0);
s=projsplx(s0);
S=diag(s);
Q=U*S*U';
q=[Q(1,1);Q(2,2);Q(1,2);q0(4:end)];
end


