% function q=Ellipse_Fitting_DR(x,nit,sigma)
%
% 
%
%
%
%
function q=Ellipse_Fitting_DR(x,nit,sigma)

%% An ellipse is parameterized as alpha x^2 + 2*beta xy + gamma y^2 + delta x + epsilon y + phi = 0
n=size(x,2);

D=zeros(6,n);
D(1,:)=x(1,:).^2;
D(2,:)=2*x(2,:).*x(1,:);
D(3,:)=x(2,:).^2;
D(4,:)=x(1,:);
D(5,:)=x(2,:);
D(6,:)=1;

K=D*D';

% %% Bias correction
% x1=x(1,:);
% x2=x(2,:);
% Delta=sigma^2*...
% [[3*n*sigma^2-6*sum(x1.^2),-6*sum(x1.*x2)                ,n*sigma^2-sum(x1.^2+x2.^2),-3*sum(x1),-sum(x2)  ,-n];...
%  [0                       ,4*n*sigma^2-4*sum(x1.^2+x2.^2),-6*sum(x1.*x2)            ,-2*sum(x2),-2*sum(x1), 0];...
%  [0                       ,0                             ,3*n*sigma^2-6*sum(x2.^2)  ,-sum(x1)  ,-3*sum(x2),-n];...
%  [0                       ,0                             ,0                         ,-n        ,0         , 0];...
%  [0                       ,0                             ,0                         ,0         ,-n        , 0];...
%  [0                       ,0                             ,0                         ,0         ,0         , 0]];
% Delta=Delta+Delta'-diag(diag(Delta));
% K=K+Delta;

%% The objective is now to solve min <q,Kq>, Tr(Q)=1, Q>=0
c1=mean(x(1,:));
c2=mean(x(2,:));
r2=(var(x(1,:))+var(x(2,:)));

u=[0.5;0;0.5;-c1;-c2;(c1^2+c2^2-r2)/2];
% And now go to the Douglas-Rachford iterative algorithm
gamma=100;  % Parameter for Douglas-Rachford in ]0,+infty[
M = gamma*K+eye(size(K));
proxf1= @(q) M\q;
proxf2= @(q) Project_on_B(q);
p=u;
for k=1:nit
    q=proxf2(p);
    p=p+1.0*(proxf1(2*q-p)-q);
end
q=proxf2(q);


function q=Project_on_B(q0)

Q0=[[q0(1),q0(2)];[q0(2),q0(3)]];
[U,S0]=eig(Q0);
s0=diag(S0);
s=projsplx(s0);
S=diag(s);
Q=U*S*U';
q=[Q(1,1);Q(2,1);Q(2,2);q0(4:end)];
end

end

