function q=Ellipse_Fitting(x,nit)

%% An ellipse is parameterized as alpha x^2 + beta xy + gamma y^2 + delta x + epsilon y + phi = 0
n=size(x,2);

D=zeros(6,n);
D(1,:)=x(1,:).^2;
D(2,:)=x(2,:).*x(1,:);
D(3,:)=x(2,:).^2;
D(4,:)=x(1,:);
D(5,:)=x(2,:);
D(6,:)=1;

K=D*D';
L=norm(K);

%% The objective is now to solve min <q,Kq>, Tr(Q)=1, Q>=0
c1=mean(x(1,:));
c2=mean(x(2,:));
r2=(var(x(1,:))+var(x(2,:)));

u=[0.5;0;0.5;-c1;-c2;c1^2+c2^2-r2];
qp=u;
q=u;

for k=1:nit
    grad=K*u;
    q=u-1/L*grad;
    
    Q0=[[q(1),q(2)];[q(2),q(3)]];
    [U,S0,V]=svd(Q0);
    s0=diag(S0);
    s=projsplx(s0);
    S=diag(s);
    Q=U*S*V';
    
    q=[Q(1);Q(2);Q(4);q(4:end)];
    u=q+0.99*(q-qp); 
    qp=q;
end

end

