addpath(genpath('/usr/share/matlab2tikz-master'))

%% Generates a set of points on an ellipse
n=150;
theta1=pi/3;

sigma=0.001;
s1=8;
s2=1;
c1=10;
c2=5;
%rng(5);

x=[];
t=linspace(-pi/16,pi/16,n);
t=[t,pi/2,pi,3*pi/2];
n=length(t);

for i=1:n
    x(1,:)=cos(theta1)*s1*cos(t)+sin(theta1)*s2*sin(t) +randn(1,n)*sigma+c1;
    x(2,:)=-sin(theta1)*s1*cos(t)+cos(theta1)*s2*sin(t) +randn(1,n)*sigma+c2;   
end

figure(1);plot(x(1,:),x(2,:),'g*');


%% Generates a fitting ellipse
s=2*max(s1,s2);
[X,Y]=meshgrid(linspace(-s+c1,s+c1,100),linspace(-s+c2,s+c2,100));

nit=50000;
for nn=nit
    tic;
    [q0,CF0]=Ellipse_Fitting_DR(x,nn);
    toc;
    tic;
    [q1,CF1]=Ellipse_Fitting_DR_SVD(x,nn);
    toc;
    tic;
    [q2,CF2]=Ellipse_Fitting_DR_EIG(x,nn);
    toc;
    
    t=linspace(0,2*pi,1000);
    figure(2);plot(x(1,:),x(2,:),'k*');axis equal
    hold on;
    title(sprintf('%i/%i',nn,nit))
    Z0=q0(1)*X.^2 + q0(2)*Y.^2 + sqrt(2)*q0(3)*X.*Y + q0(4)*X+q0(5)*Y+q0(6);
    Z1=q1(1)*X.^2 + q1(2)*Y.^2 + sqrt(2)*q1(3)*X.*Y + q1(4)*X+q1(5)*Y+q1(6);
    Z2=q2(1)*X.^2 + q2(2)*Y.^2 + sqrt(2)*q2(3)*X.*Y + q2(4)*X+q2(5)*Y+q2(6);
    contour(X,Y,Z0,[0 0],'linewidth',2,'Color',[0 0 1]);
    contour(X,Y,Z1,[0 0],'linewidth',2,'Color',[0 1 0]);
    contour(X,Y,Z2,[0 0],'linewidth',2,'Color',[1 0 0]);
    %plot(cos(theta)*s1*cos(t)+sin(theta)*s2*sin(t)+c1,-sin(theta)*s1*cos(t)+cos(theta)*s2*sin(t)+c2,'k-');
    hold off;
    legend('Data points','DR','DR-SVD','DR-EIG')
    drawnow;
end

m0=min(CF0);m1=min(CF1);m2=min(CF2);
l=length(CF1);
figure(3);semilogy(1:10:l,CF0(1:10:l)-m0,'b',1:10:l,CF1(1:10:l)-m1,'g',1:10:l,CF2(1:10:l)-m2,'r','linewidth',2);
xlabel('Iterations number');ylabel('Cost function');
legend('DR','DR-SVD','DR-EIG')

% matlab2tikz('XP_Rotation1.tex')
%  
% figure(2);
% matlab2tikz('XP_Rotation2.tex') 