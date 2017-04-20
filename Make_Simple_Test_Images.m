%% Parameters
n=128;
alpha=pi/3;
beta=pi/7;
gamma=pi/2;
s1=5;
s2=2;
s3=1;
rng(1);
center=randn(3,1)*3;
[X,Y,Z]=ndgrid(linspace(-n/2,n/2+1,n),linspace(-n/2,n/2+1,n),linspace(-n/2,n/2+1,n));


%% Points generation
A1=[[1 0 0];[0 cos(alpha) sin(alpha)];[0 -sin(alpha) cos(alpha)]];
A2=[[cos(beta) 0 sin(beta)];[0 1 0];[-sin(beta)  0 cos(beta)]];
A3=[[cos(gamma) sin(gamma) 0];[-sin(gamma) cos(gamma) 0];[0 0 1]];
R=A1*A2*A3;

for i=1:n
    XX=R(1,1)*X+R(1,2)*Y+R(1,3)*Z;
    YY=R(2,1)*X+R(2,2)*Y+R(2,3)*Z;
    ZZ=R(3,1)*X+R(3,2)*Y+R(3,3)*Z;
end

M=(XX.^2/(s1^2)+YY.^2/(s2^2)+ZZ.^2/(s3^2));
M=M/max(M(:));
A=exp(-M);

u=double(A>=0.99);
Save_3D_TIFF(uint8(u*256),'Test.tif')



% 
% sx=3;sy=10;sz=10;
% u=exp(-X.^2/(2*sx^2) - Y.^2/(2*sy^2) - Z.^2/(2*sz^2));
% u=u/max(u(:))*255;
% Save_3D_TIFF(u,'Test1.tif')
% 
% sx=10;sy=3;sz=10;
% u=exp(-X.^2/(2*sx^2) - Y.^2/(2*sy^2) - Z.^2/(2*sz^2));
% u=u/max(u(:))*255;
% Save_3D_TIFF(u,'Test2.tif')
% 
% sx=10;sy=10;sz=3;
% u=exp(-X.^2/(2*sx^2) - Y.^2/(2*sy^2) - Z.^2/(2*sz^2));
% u=u/max(u(:))*255;
% Save_3D_TIFF(u,'Test3.tif')

