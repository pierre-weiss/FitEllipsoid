%% Generates a set of points on an ellipse
n=10;
theta=0;
sigma=0.002;
s1=1;
s2=2;
c1=0;
c2=0;
rng(1);

x=[];
t=linspace(-pi/16,pi/16,n);
t=[t,pi/2,pi,3*pi/2];
n=length(t);
for i=1:n
    x(1,:)=cos(theta)*s1*cos(t)+sin(theta)*s2*sin(t) +randn(1,n)*sigma+c1;
    x(2,:)=-sin(theta)*s1*cos(t)+cos(theta)*s2*sin(t) +randn(1,n)*sigma+c2;
end

%% Generates a fitting ellipse
nit=1000;
tic;
q1=Ellipse_Fitting_DR(x,nit,sigma);
toc;
tic;
q2=Ellipse_Fitting_LLS(x);
toc;
tic;
q3=Ellipse_Fitting_ALS(x,sigma);
toc;

s=2*max(s1,s2);
[X,Y]=meshgrid(linspace(-s+c1,s+c1,100),linspace(-s+c2,s+c2,100));

figure(2);plot(x(1,:),x(2,:),'k.');axis equal

t=linspace(0,2*pi,1000);
figure(1);plot(x(1,:),x(2,:),'k.');axis equal
hold on;
Z1=q1(1)*X.^2+sqrt(2)*q1(2)*X.*Y+q1(3)*Y.^2+q1(4)*X+q1(5)*Y+q1(6);
Z2=q2(1)*X.^2+sqrt(2)*q2(2)*X.*Y+q2(3)*Y.^2+q2(4)*X+q2(5)*Y+q2(6);
Z3=q3(1)*X.^2+sqrt(2)*q3(2)*X.*Y+q3(3)*Y.^2+q3(4)*X+q3(5)*Y+q3(6);
contour(X,Y,Z1,[1e-16],'linewidth',2,'Color',[0 1 0]);
contour(X,Y,Z2,[1e-16],'linewidth',2,'Color',[1 0 0]);
contour(X,Y,Z3,[1e-16],'linewidth',2,'Color',[0 0 1]);
plot(cos(theta)*s1*cos(t)+sin(theta)*s2*sin(t)+c1,-sin(theta)*s1*cos(t)+cos(theta)*s2*sin(t)+c2,'k-');
legend('','ADMM','LLS','ALS','Ground Truth')
hold off;

%figure(2);surf(X,Y,Z);shading interp;
