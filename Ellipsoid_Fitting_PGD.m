function q=Ellipsoid_Fitting_PGD(x,nit)

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
L=norm(K);

%% The objective is now to solve min <q,Kq>, Tr(Q)=1, Q>=0
% We start with a sphere 
c1=mean(x(1,:));
c2=mean(x(2,:));
c3=mean(x(3,:));
r2=var(x(1,:))+var(x(2,:))+var(x(3,:));

u=[1/3;1/3;1/3;0;0;0;-2*c1/3;-2*c2/3;-2*c3*3;c1^2+c2^2+c3^2-r2];
qp=u;

for k=1:nit
    % Gradient descent
    grad=K*u;
    q=u-1/L*grad;
    
    % Projection
    Q0=[[q(1),q(4),q(5)];[q(4),q(2),q(6)];[q(5),q(6),q(3)]];
    [U,S0,V]=svd(Q0);
    s0=diag(S0);
    s=projsplx(s0);
    S=diag(s);
    Q=U*S*V';
    q=[Q(1,1);Q(2,2);Q(3,3);Q(2,1);Q(3,1);Q(3,2);q(7:end)];
    
    % Acceleration
    u=q+0.99*(q-qp); 
    qp=q;
end

end

