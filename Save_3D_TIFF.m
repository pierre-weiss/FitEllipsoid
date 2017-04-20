function Save_3D_TIFF(u,name)

imwrite(uint8(u(:,:,1)),name,'WriteMode','overwrite');
for i=2:size(u,3)
    imwrite(uint8(u(:,:,i)),name,'WriteMode','append');
end