addpath(genpath('/usr/share/matlab2tikz-master'))

%% Generates a set of points on an ellipse
n=10;
theta=0.1*pi;
sigma=0.05;
s1=1;
s2=2;
c1=0;
c2=0;
rng(1);

x1=[];
x2=[];
t=linspace(-pi/16,pi/16,n);
t=[t,pi/2,pi,3*pi/2];
n=length(t);

dc1=10;
for i=1:n
    x1(1,:)=cos(theta)*s1*cos(t)+sin(theta)*s2*sin(t) +randn(1,n)*sigma+c1;
    x1(2,:)=-sin(theta)*s1*cos(t)+cos(theta)*s2*sin(t) +randn(1,n)*sigma+c2;
    
    x2(1,:)=x1(1,:)+dc1;
    x2(2,:)=x1(2,:);
end

%% Generates a fitting ellipse
s=2*max(s1,s2);
[X1,Y1]=meshgrid(linspace(-s+c1,s+c1,100),linspace(-s+c2,s+c2,100));
[X2,Y2]=meshgrid(linspace(-s+c1+dc1,s+c1+dc1,100),linspace(-s+c2,s+c2,100));

nit=50000;
for nn=1:10000:nit
    tic;
    [q1,CF1]=Ellipse_Fitting_DR(x1,nn,sigma);
    toc;
    tic;
    [q2,CF2]=Ellipse_Fitting_DR(x2,nn,sigma);
    toc;
    
    t=linspace(0,2*pi,1000);
    figure(2);plot(x1(1,:),x1(2,:),'k*');axis equal
    title(sprintf('%i/%i',nn,nit))
    hold on;
    Z1=q1(1)*X1.^2 + q1(2)*Y1.^2 + sqrt(2)*q1(3)*X1.*Y1 + q1(4)*X1+q1(5)*Y1+q1(6);
    Z2=q2(1)*X2.^2 + q2(2)*Y2.^2 + sqrt(2)*q2(3)*X2.*Y2 + q2(4)*X2+q2(5)*Y2+q2(6);
    contour(X1,Y1,Z1,[0 0],'linewidth',2,'Color',[0 1 0]);
    contour(X2-dc1,Y2,Z2,[0 0],'linewidth',2,'Color',[1 0 0]);
    plot(cos(theta)*s1*cos(t)+sin(theta)*s2*sin(t)+c1,-sin(theta)*s1*cos(t)+cos(theta)*s2*sin(t)+c2,'k-');
    hold off;
    legend('Data points','ADMM1','ADMM2','Ground Truth')
    drawnow;
end

m=min(CF1);l=length(CF1);figure(1);semilogy(1:10:l,CF1(1:10:l)-m,'r',1:10:l,CF2(1:10:l)-m,'b','linewidth',2);xlabel('Iterations number');ylabel('Cost function');legend('Centered','Shifted')
% matlab2tikz('Compare_CostFunction.tex')
% 
% figure(2);
% matlab2tikz('XP_Shift_100it.tex')
% 