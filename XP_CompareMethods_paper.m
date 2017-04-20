addpath(genpath('/usr/share/matlab2tikz-master'))

%% Generates a set of points on an ellipse
n=15;
theta1=pi/3;

%sigma=[0,0.1,0.3];
sigma=[0.1,0.15,0.2,0.25,0.3,0.35];
s1=4;
s2=1;
c1=4;
c2=5;

x=[];

%% Definition of sampling schemes
t=cell(4,1);
t{1}=rand(1,n)*2*pi;

t{2}=linspace(-pi/16,pi/16,n);
t{2}=[t{2},pi/2,pi,3*pi/2];

t{3}=linspace(-pi/16,pi/16,n);
t{3}=[t{3},pi/2,pi];

t{4}=linspace(-pi/4+pi/2,pi/4+pi/2,n);


for j=4:4
    n=length(t{j});
    for k=1:length(sigma)
        %% Generates the input coordinates
        x=zeros(2,n);
        x(1,:)=cos(theta1)*s1*cos(t{j})+sin(theta1)*s2*sin(t{j}) +randn(1,n)*sigma(k)+c1;
        x(2,:)=-sin(theta1)*s1*cos(t{j})+cos(theta1)*s2*sin(t{j}) +randn(1,n)*sigma(k)+c2;
        
        %% Generates a fitting ellipse
        s=1.2*max(s1,s2);
        Rx=[-s+c1,s+c1];
        Ry=[-s+c2,s+c2];
        
        nit=50000;
        for nn=nit
            tic;
            [q0,CF0]=Ellipse_Fitting_DR(x,nn);
            toc;
            tic;
            [q1,CF1]=Ellipse_Fitting_DR_SVD(x,nn);
            toc;
            tic;
            q2=Ellipse_Fitting_LLS(x);
            toc;
            tic;
            q3=Ellipse_Fitting_LLS_SVD(x);
            toc;
            
            tt=linspace(0,2*pi,1000);
            figure(2);hold off;
            plot(x(1,:),x(2,:),'k*');axis equal
            hold on;
            title(sprintf('Iteration %i/%i',nn,nit))
            
            DisplayEllipse(Rx,Ry,q0,[0 0 1]);
            DisplayEllipse(Rx,Ry,q1,[0 1 0]);
            DisplayEllipse(Rx,Ry,q2,[1 0 0]);
            DisplayEllipse(Rx,Ry,q3,[1 1 0]);
            
            plot(cos(theta1)*s1*cos(tt)+sin(theta1)*s2*sin(tt)+c1,-sin(theta1)*s1*cos(tt)+cos(theta1)*s2*sin(tt)+c2,'k--');
            hold off;
            legend('Data points','DR','DR-SVD','LLS','LLS-SVD','Ground truth');
            drawnow;
        end
        
        m0=min(CF0);m1=min(CF1);
        l=50000;
        figure(3);semilogy(1:l,CF0(1:l)-m0,'b',1:l,CF1(1:l)-m1,'g');
        xlabel('Iterations number');ylabel('Cost function');
        legend('DR','DR-SVD')
        
%         matlab2tikz(sprintf('XP_Compare1_%i_%i.tex',j,k));
%         figure(2);
%         matlab2tikz(sprintf('XP_Compare2_%i_%i.tex',j,k));
        pause;
    end
end

