% function q=Ellipse_Fitting_ALS(x,sigma)
% 
% This function 
%
%
%
function q=Ellipse_Fitting_ALS(x,sigma)

%% An ellipse is parameterized as alpha x^2 + 2 * beta xy + gamma y^2 + delta x + epsilon y + phi = 0
n=size(x,2);

D=zeros(6,n);
D(1,:)=x(1,:).^2;
D(2,:)=2*x(2,:).*x(1,:);
D(3,:)=x(2,:).^2;
D(4,:)=x(1,:);
D(5,:)=x(2,:);
D(6,:)=1;

K=D*D';
x1=x(1,:);
x2=x(2,:);

Delta=sigma^2*...
[[3*n*sigma^2-6*sum(x1.^2),-6*sum(x1.*x2)                ,n*sigma^2-sum(x1.^2+x2.^2),-3*sum(x1),-sum(x2)  ,-n];...
 [0                       ,4*n*sigma^2-4*sum(x1.^2+x2.^2),-6*sum(x1.*x2)            ,-2*sum(x2),-2*sum(x1), 0];...
 [0                       ,0                             ,3*n*sigma^2-6*sum(x2.^2)  ,-sum(x1)  ,-3*sum(x2),-n];...
 [0                       ,0                             ,0                         ,-n        ,0         , 0];...
 [0                       ,0                             ,0                         ,0         ,-n        , 0];...
 [0                       ,0                             ,0                         ,0         ,0         , 0]];
Delta=Delta+Delta'-diag(diag(Delta));

[V,D]=svd(K+Delta);
q=V(:,end);

disp(diag(D))
