clear all;
xx=dlmread('/media/weiss/Donnees/Works/Workspace_Eclipse/FitEllipsoid/dataPoints.txt');
%xr=dlmread('/media/weiss/Donnees/Works/Workspace_Eclipse/FitEllipsoid/Points_Invariant.txt');

n=size(xx,2);
t=mean(xx,2);
xb(1,:)=xx(1,:)-t(1);
xb(2,:)=xx(2,:)-t(2);
xb(3,:)=xx(3,:)-t(3);
[U,S]=eig(xb*xb');
sp=max(diag(S),1e-15*ones(3,1));
P=diag(sp.^(-0.5))*U';

x=P*xb;

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
nit=113;
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

%% Go back to the initial basis. The ellipsoid is given by <Ax,x>+<b,x>+c=0
A=P'*A2*P;
b=-2*A*t+P'*b2;
c=(A2*P*t)'*(P*t)-b2'*P*t+c2;

n=trace(A);
A=A/n;
b=b/n;
c=c/n;

q=[A(1,1);A(2,2);A(3,3);A(2,1);A(3,1);A(3,2);b(1);b(2);b(3);c];
%q=q/(A(1,1)+A(2,2)+A(3,3)); % Just a normalization to stay on the simplex


z=A\(-b/2);
[U,S,V]=svd(A);
r2=z'*A*z-c;
l=sqrt(r2./diag(S));


