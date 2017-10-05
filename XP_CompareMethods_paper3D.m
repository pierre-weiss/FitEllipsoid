addpath(genpath('/usr/share/matlab2tikz-master'))

%% Generates a set of points on an ellipse
n=50;
alpha=pi/3;
beta=-pi/7;
gamma=pi/2;

s1=5;
s2=2;
s3=3;

smax=max([s1,s2,s3]);

rng(1);
center=0*[-19;35;22];
% options
affiche=1;
partiel=1;

% distribution partielle des points avec un cosinus avec la direction
% donnee qui est plus grand que t(j)
t=[0.9 0.8 0 -2];
direction=[1,-2,3];
direction=direction/norm(direction);

sigma=[0 0.01 0.1];
 
for j=3%:4
    
    for k=1:length(sigma)
        %% Points generation
        A1=[[1 0 0];[0 cos(alpha) sin(alpha)];[0 -sin(alpha) cos(alpha)]];
        A2=[[cos(beta) 0 sin(beta)];[0 1 0];[-sin(beta)  0 cos(beta)]];
        A3=[[cos(gamma) sin(gamma) 0];[-sin(gamma) cos(gamma) 0];[0 0 1]];
        R=A1*A2*A3;
        
        
        
        xx=zeros(3,n);
        kk=1;
        while kk<=n
            v=randn(3,1);
            v=v/norm(v);
            if dot(v,direction)>t(j)
                xx(:,kk)=v;
                kk=kk+1;
            end
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
        %%
        
        x=x+randn(size(x))*sigma(k)+repmat(center,1,n);
        
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
            q3=Ellipsoid3D_Fitting_LLS(x);
            toc;
            tic;
            q4=Ellipsoid3D_Fitting_LLS_SVD(x);
            toc;
            
            r=smax;
            Rx=[center(1)-smax-0.1*r center(1)+smax+0.1*r];
            Ry=[center(2)-smax-0.1*r center(2)+smax+0.1*r];
            Rz=[center(3)-smax-0.1*r center(3)+smax+0.1*r];
            
            no=1+10*j+100*k;
            figure(no);plot3(x(1,:),x(2,:),x(3,:),'.');
            axis equal
            hold on;
            DisplayEllipsoid(Rx,Ry,Rz,q0,[1 0 0]);
            
            figure(no+1);plot3(x(1,:),x(2,:),x(3,:),'.');
            axis equal
            hold on;
            DisplayEllipsoid(Rx,Ry,Rz,q1,[0 1 0]);
            
            figure(no+2);plot3(x(1,:),x(2,:),x(3,:),'.');
            axis equal
            hold on;
            DisplayEllipsoid(Rx,Ry,Rz,q2,[0 0 1]);
            
            figure(no+3);plot3(x(1,:),x(2,:),x(3,:),'.');
            axis equal
            hold on;
            DisplayEllipsoid(Rx,Ry,Rz,q3,[1 1 0]);
            
            figure(no+4);plot3(x(1,:),x(2,:),x(3,:),'.');
            axis equal
            hold on;
            DisplayEllipsoid(Rx,Ry,Rz,q4,[1 0 1]);
            
        end
        
        m0=min(CF0);m1=min(CF1);m2=min(CF2);
        l=nit;
        figure(no+1);semilogy(1:l,CF0(1:l)-m0,'r',1:l,CF1(1:l)-m1,'g',1:l,CF2(1:l)-m2,'b');
        xlabel('Iterations number');ylabel('Cost function');
        legend('DR','DR-SVD','oblique')
        
        %         matlab2tikz(sprintf('XP_Compare1_%i_%i.tex',j,k));
        %         figure(2);
        %         matlab2tikz(sprintf('XP_Compare2_%i_%i.tex',j,k));
        pause;
    end
end

