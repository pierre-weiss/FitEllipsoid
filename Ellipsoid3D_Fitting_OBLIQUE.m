function [qs,CF,A,b,c]=Ellipsoid3D_Fitting_OBLIQUE(x,nit)

%% An ellipse is parameterized as
%% a1 x^2 + a2 y^2 + a3 z^2 + a4 xy + a5 xz + a6 yz + a7 x + a8 y + a9 z + a10 = 0
% in this script we minimize the cost function <q,Kq>, s.t. Tr(Q)=1, Q>=0
%
% we use a preconditioning that consists of changing the basis
% for the data - the initial guess is the smallest eigenvector of K

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
[~,s1,v1]=svd(K);
a=v1(:,10);
a=a/sum(a(1:3));
if a(1)<0
    a=-a;
end

A=[a(1) a(4)/sqrt(2) a(5)/sqrt(2);a(4)/sqrt(2) a(2) a(6)/sqrt(2);a(5)/sqrt(2) a(6)/sqrt(2) a(3)];
if min(eig(A))<0
    disp('Beware, the OBLIQUE solution is not an ellipse');
end

sm=[a(7);a(8);a(9)];
c=-0.5*(A\sm);
r=sqrt(1e-16-a(10)+dot(A*c,c));
% dans la nouvelle base on a xx
[uA,sA,vA]=svd(A);
Ademi=uA*(max(sA,1e-16).^0.5)*vA';
xx=1/r*Ademi*(x-repmat(c,1,n));

D=zeros(10,n);
D(1,:)=xx(1,:).^2;
D(2,:)=xx(2,:).^2;
D(3,:)=xx(3,:).^2;
D(4,:)=sqrt(2)*xx(1,:).*xx(2,:);
D(5,:)=sqrt(2)*xx(1,:).*xx(3,:);
D(6,:)=sqrt(2)*xx(2,:).*xx(3,:);
D(7,:)=xx(1,:);
D(8,:)=xx(2,:);
D(9,:)=xx(3,:);
D(10,:)=1;

K=D*D';
%% The objective is now to solve min <q,Kq>, Tr(Q)=1, Q>=0 using a Douglas-Rachford algorithm
% We start with a sphere
c1=mean(xx(1,:));
c2=mean(xx(2,:));
c3=mean(xx(3,:));
r2=var(xx(1,:))+var(xx(2,:))+var(xx(3,:));
u=[1/3;1/3;1/3;0;0;0;-2*c1/3;-2*c2/3;-2*c3/3;(c1^2+c2^2+c3^2-r2)/3];

% And now go to the Douglas-Rachford iterative algorithm
gamma=10; % Parameter for Douglas-Rachford in ]0,+infty[
M = gamma*K+eye(size(K));
proxf1= @(q) M\q;
proxf2= @(q) Project_on_B(q);
[u1,s1,v1]=svd(K);
p=real(v1(:,10));
p=p/sum(p(1:3));



CF=zeros(nit+1,1);

for k=1:nit
    
    q=proxf2(p);
    p=p+proxf1(2*q-p)-q;
    CF(k)=0.5*q'*K*q;
end
q=proxf2(q);
q2=q;
CF(nit+1)=0.5*q'*K*q;


qs=zeros(10,1);
% back to the original coordinates
B=[q2(1) q2(4)/sqrt(2) q2(5)/sqrt(2);q2(4)/sqrt(2) q2(2) q2(6)/sqrt(2);...
    q2(5)/sqrt(2) q2(6)/sqrt(2) q2(3)];
sm=[q2(7);q2(8);q2(9)];
cprime=-0.5*(B\sm);
rprime=sqrt(-q2(10)+dot(B*cprime,cprime));
qs=zeros(10,1);
M=Ademi*B*Ademi;
qs(1:6)=1/(r^2)*[M(1,1);M(2,2);M(3,3);sqrt(2)*M(2,1);sqrt(2)*M(3,1);sqrt(2)*M(3,2)];
qs(7:9)=-(2/r)*Ademi*B*(1/r*Ademi*c+cprime);
qs(10)=dot(B*(1/r*Ademi*c+cprime),1/r*Ademi*c+cprime)-rprime^2;
qs=qs/sum(qs(1:3));
qs=real(qs);

end

function q=Project_on_B(q0)
% ici on remplace la projection sur B par un rescaling qui fait intervenir
% les coefficients diagonaux ainsi que le rayon
Q0=[[q0(1),q0(4)/sqrt(2),q0(5)/sqrt(2)];[q0(4)/sqrt(2),q0(2),q0(6)/sqrt(2)];...
    [q0(5)/sqrt(2),q0(6)/sqrt(2),q0(3)]];
[U,S0,V]=svd(Q0);
sm=[q0(7);q0(8);q0(9)];
c=-0.5*(Q0\sm);
r0=sqrt(-q0(10)+dot(Q0*c,c));

s0=diag(S0);
if min(s0)<=0
    disp 'tata'
    s=projsplx(s0);
    S=diag(s);
    Q=U*S*U';
    q=zeros(10,1);
    q(1:6)=[Q(1,1);Q(2,2);Q(3,3);sqrt(2)*Q(2,1);sqrt(2)*Q(3,1);sqrt(2)*Q(3,2)];
    q(7:9)=-2*Q*c;
    q(10)=dot(Q*c,c)-r0^2;
else
    e1=[1;1;1;0;0;0;0;0;0;-1];
    e1=e1/norm(e1);
    q=q0-(dot(q0,e1)-1)*e1;
end
end