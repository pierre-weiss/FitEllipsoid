addpath(genpath('/usr/share/matlab2tikz-master'))

n=4;
%rng(5);
%x=rand(2,n);

x=[[0.6;0.1],[0.2;0.5],[1;0.5],[0.5;1]];

%
D=zeros(6,n);
D(1,:)=x(1,:).^2;
D(2,:)=x(2,:).^2;
D(3,:)=sqrt(2)*x(2,:).*x(1,:);
D(4,:)=x(1,:);
D(5,:)=x(2,:);
D(6,:)=1;
K=D*D';

[~,S,V1]=svd(K);
q1=V1(:,6);
q1=q1/sum(q1(1:2));

q2=V1(:,5);
q2=q2/sum(q2(1:2));


figure(1);hold off;
for alpha=0.1:0.1:0.7

    q=alpha*q1+(1-alpha)*q2;

    A=[[q(1),q(3)/sqrt(2)];[q(3)/sqrt(2),q(2)]];
    b=[q(4);q(5)];
    c=q(6);
    
    DisplayEllipse([-1 2],[-1 2],q,[alpha/0.8 1-(alpha/0.8) 0]);
    hold on;
    plot(x(1,:),x(2,:),'b*');
    pause(0.1);
end
axis equal;

matlab2tikz('XP_NonUniqueness.tex')
