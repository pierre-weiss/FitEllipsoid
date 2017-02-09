% function q=Ellipse_Fitting_DR_EIG(x,nit,sigma)
%
function [q,CF]=Ellipse_Fitting_DR_EIG(x,nit)

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

%% Change of basis
[~,~,V1]=svd(K);
a=V1(:,6);
a=a/sum(a(1:2));
if a(1)<0
    a=-a;
end

A=[[a(1) a(3)/sqrt(2)];[a(3)/sqrt(2),a(2)]];
b=[a(4);a(5)];
c=-0.5*(A\b);
r=sqrt(-a(6)+dot(A*c,c));
% dans la nouvelle base on a xx
[uA,sA,vA]=svd(A);
Ademi=uA*(max(sA,0).^0.5)*vA';

P=1/r*Ademi;
t=c;

x=P*(x-repmat(c,1,n));


%% Redefines K with the new basis
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
gamma=1;  % Parameter for Douglas-Rachford in ]0,+infty[
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

    function q=Project_on_B(q0)
        Q0=[[q0(1),q0(3)/sqrt(2)];[q0(3)/sqrt(2),q0(2)]];
        [U,S0]=eig(Q0);
        s0=diag(S0);
        s=projsplx(s0);
        S=diag(s);
        Q=U*S*U';
        q=[Q(1,1);Q(2,2);sqrt(2)*Q(1,2);q0(4:end)];
    end

end

