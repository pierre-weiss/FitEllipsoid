%% Parameters
nit=1000;
n=1000;
alpha=0;%pi/3;
beta=0;%pi/7;
gamma=0;%pi/2;
c1=2;
c2=-1;
s1=1;
s2=2;
s3=10;
sigma=0; %noise variance
rng(1);
center=10*randn(3,1);


%% Points generation
A1=[[1 0 0];[0 cos(alpha) sin(alpha)];[0 -sin(alpha) cos(alpha)]];
A2=[[cos(beta) 0 sin(beta)];[0 1 0];[-sin(beta)  0 cos(beta)]];
A3=[[cos(gamma) sin(gamma) 0];[-sin(gamma) cos(gamma) 0];[0 0 1]];
R=A1*A2*A3;
xx=randn(3,n);
for i=1:n
    xx(:,i)=xx(:,i)/norm(xx(:,i));
end
xx(1,:)=xx(1,:)*s1;
xx(2,:)=xx(2,:)*s2;
xx(3,:)=xx(3,:)*s3;
x=zeros(size(xx));
for i=1:n
    x(1,:)=R(1,1)*xx(1,:)+R(1,2)*xx(2,:)+R(1,3)*xx(3,:);
    x(2,:)=R(2,1)*xx(1,:)+R(2,2)*xx(2,:)+R(2,3)*xx(3,:);
    x(3,:)=R(3,1)*xx(1,:)+R(3,2)*xx(2,:)+R(3,3)*xx(3,:);
end
x=x+randn(size(x))*sigma+repmat(center,1,n);

%% Generates a fitting ellipse
tic;
q=Ellipsoid_Fitting_Centering(x,nit);
toc;

%% Displays the result
s=1.5*max(abs(x(:)));
[X,Y,Z]=meshgrid(linspace(-s,s,100),linspace(-s,s,100),linspace(-s,s,100));

figure(1);plot3(x(1,:),x(2,:),x(3,:),'.');axis equal
hold on;
V=q(1)*X.^2+q(2)*Y.^2+q(3)*Z.^2+q(4)*X.*Y+q(5)*X.*Z+q(6)*Y.*Z+q(7)*X+q(8)*Y+q(9)*Z+q(10);
p=patch(isosurface(X,Y,Z,V,0));
isonormals(X,Y,Z,V,p);
set(p,'FaceColor','red','EdgeColor','none');
set(p,'FaceAlpha',0.5);
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud
hold off;

%% Finding the matrices
A0=[[q(1) q(4)/2 q(5)/2];[q(4)/2 q(2) q(6)/2];[q(5)/2 q(6)/2 q(3)]];
c=-A0\[q(7)/2;q(8)/2;q(9)/2];

alpha= 1/(c'*(A0*c) - q(10));


% prendre racine réelle alpha puis le code continue comme suit:
A=A0*alpha;
[U,S]=eig(A);
disp(1./sqrt(diag(S)))





