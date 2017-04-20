function q=Ellipsoid_Fitting(x,nit)

%% An ellipse is parameterized as
%% a1 x^2 + a2 y^2 + a3 z^2 + a4 xy + a5 xz + a6 yz + a7 x + a8 y + a9 z + a10 = 0
n=size(x,2);

D=zeros(10,n);
D(1,:)=x(1,:).^2;
D(2,:)=x(2,:).^2;
D(3,:)=x(3,:).^2;
D(4,:)=x(1,:).*x(2,:);
D(5,:)=x(1,:).*x(3,:);
D(6,:)=x(2,:).*x(3,:);
D(7,:)=x(1,:);
D(8,:)=x(2,:);
D(9,:)=x(3,:);
D(10,:)=1;

K=D*D';

%% The objective is now to solve min <q,Kq>, Tr(Q)=1, Q>=0 using a Douglas-Rachford algorithm
% We start with a sphere
c1=mean(x(1,:));
c2=mean(x(2,:));
c3=mean(x(3,:));
r2=var(x(1,:))+var(x(2,:))+var(x(3,:));
u=[1/3;1/3;1/3;0;0;0;-2*c1/3;-2*c2/3;-2*c3/3;(c1^2+c2^2+c3^2-r2)/3];

% And now go to the Douglas-Rachford iterative algorithm
gamma=0.01; % Parameter for Douglas-Rachford in ]0,+infty[
% Inverse = inv(gamma*K+eye(size(K)));
% proxf1= @(q) Inverse*q;
M = gamma*K+eye(size(K));
proxf1= @(q) M\q;
proxf2= @(q) Project_on_B(q);
p=u;
for k=1:nit
    q=proxf2(p);
    p=p+proxf1(2*q-p)-q;
end
q=proxf2(q);

end

function q=Project_on_B(q0)

Q0=[[q0(1),q0(4)/2,q0(5)/2];[q0(4)/2,q0(2),q0(6)/2];[q0(5)/2,q0(6)/2,q0(3)]];
[U,S0]=eig(Q0);
s0=diag(S0);
s=projsplx(s0);
S=diag(s);
Q=U*S*U';
q=[Q(1,1);Q(2,2);Q(3,3);2*Q(2,1);2*Q(3,1);2*Q(3,2);q0(7:end)];
end