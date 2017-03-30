function q=Ellipse_Fitting_LLS(x)

%% An ellipse is parameterized as alpha x^2 + beta xy + gamma y^2 + delta x + epsilon y + phi = 0
n=size(x,2);

D=zeros(6,n);
D(1,:)=x(1,:).^2;
D(2,:)=2*x(2,:).*x(1,:);
D(3,:)=x(2,:).^2;
D(4,:)=x(1,:);
D(5,:)=x(2,:);
D(6,:)=1;

K=D*D';

[V,D]=svd(K);
q=V(:,end);

disp(diag(D))
