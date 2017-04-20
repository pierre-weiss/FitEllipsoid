addpath(genpath('/usr/share/matlab2tikz-master'))

%% Generates a set of points on an ellipse
n=20;
alpha=pi/3;
beta=-pi/7;
gamma=pi/2;

s1=5;
s2=2;
s3=3;

rng(22);
center=[10;0;0];
% options
affiche=1;
partiel=1;

% distribution partielle des points avec un cosinus avec la direction
% donnee qui est plus grand que t(j)
%t=[0.9 0.8 0 -2];
%direction=[1,-2,3];
%direction=direction/norm(direction);

%sigma=[0 0.01 0.1];
sigma=[0.2];



for k=1:length(sigma)
    %% Points generation
    A1=[[1 0 0];[0 cos(alpha) sin(alpha)];[0 -sin(alpha) cos(alpha)]];
    A2=[[cos(beta) 0 sin(beta)];[0 1 0];[-sin(beta)  0 cos(beta)]];
    A3=[[cos(gamma) sin(gamma) 0];[-sin(gamma) cos(gamma) 0];[0 0 1]];
    R=A1*A2*A3;
    Sigma=diag([s1,s2,s3]);
    
    %% Draw points at random on equators
    xx=zeros(3,n);
    kk=1;
    while kk<=n
        v=randn(2,1);
        if mod(kk,3)==0
            v=[v;0]/norm(v);
        end
        if mod(kk,3)==1
            v=[v(1);0;v(2)]/norm(v);
        end
        if mod(kk,3)==2
            v=[0;v]/norm(v);
        end
        w=Sigma\(R\v);
        xx(:,kk)=w/norm(w);
        kk=kk+1;
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
    
    %% Compute equators
    nbe=100;
    xe1=zeros(3,nbe);
    xe2=zeros(3,nbe);
    xe3=zeros(3,nbe);
    ind=1;
    for kk=0:3*nbe-1
        t=2*pi/(3*nbe)*kk;
        v=randn(2,1);
        if mod(kk,3)==0
            xe1(1,ind)=0*s1;
            xe1(2,ind)=cos(t);
            xe1(3,ind)=sin(t);
            xe1(:,ind)=Sigma\(R\xe1(:,ind));
            xe1(:,ind)=xe1(:,ind)/norm(xe1(:,ind));
        end
        if mod(kk,3)==1
            xe2(1,ind)=cos(t);
            xe2(2,ind)=0;
            xe2(3,ind)=sin(t);
            xe2(:,ind)=Sigma\(R\xe2(:,ind));
            xe2(:,ind)=xe2(:,ind)/norm(xe2(:,ind));
        end
        if mod(kk,3)==2
            xe3(1,ind)=cos(t);
            xe3(2,ind)=sin(t);
            xe3(3,ind)=0;
            xe3(:,ind)=Sigma\(R\xe3(:,ind));
            xe3(:,ind)=xe3(:,ind)/norm(xe3(:,ind));

            ind=ind+1;
        end
    end
    xe1(1,:)=xe1(1,:)*s1;xe1(2,:)=xe1(2,:)*s2;xe1(3,:)=xe1(3,:)*s3;
    xe2(1,:)=xe2(1,:)*s1;xe2(2,:)=xe2(2,:)*s2;xe2(3,:)=xe2(3,:)*s3;
    xe3(1,:)=xe3(1,:)*s1;xe3(2,:)=xe3(2,:)*s2;xe3(3,:)=xe3(3,:)*s3;
    x1=zeros(size(xe1));
    x2=zeros(size(xe1));
    x3=zeros(size(xe1));
    for i=1:n
        x1(1,:)=R(1,1)*xe1(1,:)+R(1,2)*xe1(2,:)+R(1,3)*xe1(3,:);
        x1(2,:)=R(2,1)*xe1(1,:)+R(2,2)*xe1(2,:)+R(2,3)*xe1(3,:);
        x1(3,:)=R(3,1)*xe1(1,:)+R(3,2)*xe1(2,:)+R(3,3)*xe1(3,:);
        
        x2(1,:)=R(1,1)*xe2(1,:)+R(1,2)*xe2(2,:)+R(1,3)*xe2(3,:);
        x2(2,:)=R(2,1)*xe2(1,:)+R(2,2)*xe2(2,:)+R(2,3)*xe2(3,:);
        x2(3,:)=R(3,1)*xe2(1,:)+R(3,2)*xe2(2,:)+R(3,3)*xe2(3,:);

        x3(1,:)=R(1,1)*xe3(1,:)+R(1,2)*xe3(2,:)+R(1,3)*xe3(3,:);
        x3(2,:)=R(2,1)*xe3(1,:)+R(2,2)*xe3(2,:)+R(2,3)*xe3(3,:);
        x3(3,:)=R(3,1)*xe3(1,:)+R(3,2)*xe3(2,:)+R(3,3)*xe3(3,:);        
    end
    
    
    %% Adds noise
    kk=1;
    while kk<=n
        v=randn(2,1);
        if mod(kk,3)==0
            x(1:2,kk)=x(1:2,kk) + v*sigma(k);
        end
        if mod(kk,3)==1
            x([1 3],kk)=x([1 3],kk) + v*sigma(k);
        end
        if mod(kk,3)==2
            x(2:3,kk)=x(2:3,kk) + v*sigma(k);
        end
        kk=kk+1;
    end

    
    nit=1e3;
    for nn=nit
        tic;
        [q0,CF0]=Ellipsoid3D_Fitting_DR(x,nn);
        toc;
        tic;
        [q1,CF1]=Ellipsoid3D_Fitting_DR_SVD(x,nn);
        toc;
        tic;
        [q2,CF2]=Ellipsoid3D_Fitting_OBLIQUE(x,nn);
        toc;
        tic;
        q3=Ellipsoid_Fitting_LLS(x);
        toc;
        tic;
        q4=Ellipsoid3D_Fitting_LLS_SVD(x);
        toc;
        
        r=max([s1 s2 s3]);
        Rx=[-1.3*r 1.3*r];
        Ry=[-1.3*r 1.3*r];
        Rz=[-1.3*r 1.3*r];
        
        no=1+100*k;
        figure(no);plot3(x(1,:),x(2,:),x(3,:),'.k');
        axis equal
        hold on;
        plot3(x1(1,:),x1(2,:),x1(3,:),'-b');
        plot3(x2(1,:),x2(2,:),x2(3,:),'-b');
        plot3(x3(1,:),x3(2,:),x3(3,:),'-b');
        DisplayEllipsoid(Rx,Ry,Rz,q0,[1 0 0]);
        title('DR')
        
        fig=figure(no+1);plot3(x(1,:),x(2,:),x(3,:),'.k');
        axis equal
        hold on;
        plot3(x1(1,:),x1(2,:),x1(3,:),'-b');
        plot3(x2(1,:),x2(2,:),x2(3,:),'-b');
        plot3(x3(1,:),x3(2,:),x3(3,:),'-b');        
        DisplayEllipsoid(Rx,Ry,Rz,q1,[0 1 0]);
        %title('DR-SVD')

        
        figure(no+2);plot3(x(1,:),x(2,:),x(3,:),'.k');
        axis equal
        hold on;
        plot3(x1(1,:),x1(2,:),x1(3,:),'-b');
        plot3(x2(1,:),x2(2,:),x2(3,:),'-b');
        plot3(x3(1,:),x3(2,:),x3(3,:),'-b');
        DisplayEllipsoid(Rx,Ry,Rz,q2,[0 0 1]);
        title('Oblique')

        
        figure(no+3);plot3(x(1,:),x(2,:),x(3,:),'.k');
        axis equal
        hold on;
        plot3(x1(1,:),x1(2,:),x1(3,:),'-b');
        plot3(x2(1,:),x2(2,:),x2(3,:),'-b');
        plot3(x3(1,:),x3(2,:),x3(3,:),'-b');
        DisplayEllipsoid(Rx,Ry,Rz,q3,[1 1 0]);
        title('LLS')

        
        figure(no+4);plot3(x(1,:),x(2,:),x(3,:),'.k');
        axis equal
        hold on;
        plot3(x1(1,:),x1(2,:),x1(3,:),'-b');
        plot3(x2(1,:),x2(2,:),x2(3,:),'-b');
        plot3(x3(1,:),x3(2,:),x3(3,:),'-b');
        DisplayEllipsoid(Rx,Ry,Rz,q4,[1 0 1]);
        title('LLS-SVD')
    end
    
    m0=min(CF0);m1=min(CF1);m2=min(CF2);
    l=nit;
    figure(no+10);semilogy(1:l,CF0(1:l)-m0,'r',1:l,CF1(1:l)-m1,'g',1:l,CF2(1:l)-m2,'b');
    xlabel('Iterations number');ylabel('Cost function');
    legend('DR','DR-SVD','oblique')
        
    %matlab2tikz('Ellipsoid3D_1.tex');
end


figure(fig);
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,sprintf('XP_3D_%1.2e_%i.png',sigma,n),'-dpng');



