addpath(genpath('/usr/share/matlab2tikz-master'))

%% Generates a set of points on an ellipse
n=10;
theta1=0;

L=[[8,0];[0,1]];

sigma=0.0;
s1=2;
s2=1;
c1=0;
c2=0;
rng(1);

x1=[];
x2=[];
t=linspace(-pi/16,pi/16,n);
t=[t,pi/2,pi,3*pi/2];
n=length(t);

for i=1:n
    x1(1,:)=cos(theta1)*s1*cos(t)+sin(theta1)*s2*sin(t) +randn(1,n)*sigma+c1;
    x1(2,:)=-sin(theta1)*s1*cos(t)+cos(theta1)*s2*sin(t) +randn(1,n)*sigma+c2;
    
    x2=L*x1+[5;10];
end

figure(1);plot(x1(1,:),x1(2,:),'g*');
hold on;
plot(x2(1,:),x2(2,:),'r*');axis equal;
legend('x1','x2')
hold off;


%% Generates a fitting ellipse
s=2*max(s1*d1,s2*d2);
[X1,Y1]=meshgrid(linspace(-s+c1,s+c1,100),linspace(-s+c2,s+c2,100));
X2=X1;
Y2=Y1;
%X2=cos(theta2)*d1*X1+sin(theta2)*d2*Y1;
%Y2=-sin(theta2)*d1*X1+cos(theta2)*d2*Y1;

nit=1000;
for nn=nit
    tic;
    [q1,CF1]=Ellipse_Fitting_DR_SVD(x1,nn);
    toc;
    tic;
    [q2,CF2]=Ellipse_Fitting_DR_SVD(x2,nn);
    toc;
    
    t=linspace(0,2*pi,1000);
    figure(2);plot(x1(1,:),x1(2,:),'k*');axis equal
    hold on;
    figure(2);plot(x2(1,:),x2(2,:),'ko');axis equal
    title(sprintf('%i/%i',nn,nit))
    Z1=q1(1)*X1.^2 + q1(2)*Y1.^2 + sqrt(2)*q1(3)*X1.*Y1 + q1(4)*X1+q1(5)*Y1+q1(6);
    Z2=q2(1)*X2.^2 + q2(2)*Y2.^2 + sqrt(2)*q2(3)*X2.*Y2 + q2(4)*X2+q2(5)*Y2+q2(6);
    contour(X1,Y1,Z1,[0 0],'linewidth',2,'Color',[0 1 0]);
    contour(X2,Y2,Z2,[0 0],'linewidth',2,'Color',[1 0 0]);
    %plot(cos(theta)*s1*cos(t)+sin(theta)*s2*sin(t)+c1,-sin(theta)*s1*cos(t)+cos(theta)*s2*sin(t)+c2,'k-');
    hold off;
    legend('Data points 1','Data points 2','ADMM1','ADMM2')
    drawnow;
end

m1=min(CF1);m2=min(CF2);l=length(CF1);figure(3);semilogy(1:10:l,CF1(1:10:l)-m1,'r',1:10:l,CF2(1:10:l)-m2,'b','linewidth',2);
xlabel('Iterations number');ylabel('Cost function');
legend('Circle','Ellipse')

% matlab2tikz('XP_Rotation1.tex')
%  
% figure(2);
% matlab2tikz('XP_Rotation2.tex') 