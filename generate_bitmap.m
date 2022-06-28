function [image_bitmap] = generate_bitmap(coords,xlength,ylength,xdim,ydim)

image = 'fix_this_test2';

x0 = coords(1);
y0 = coords(2);

%Read image and convert to binary
image = imread(strcat(image,'.png'));
image = im2bw(image);

figure(),
imshow(image);
title('Black and White Bitmap, pre-conversion');

%Resize and create array
image_reshape = imresize(image,[xdim,ydim]);

xside=linspace(x0,x0+xlength,xdim);
yside=linspace(y0,y0+ylength,ydim);
bitmap = zeros(size(image_reshape));

for ii=1:xdim
    x=xside(ii);
    for jj=1:ydim
        
        y=yside(jj);
        pixel=image_reshape(ii,jj);
        
        if pixel==0
            bitmap(ii,jj) = 1;
        elseif pixel==1
            bitmap(ii,jj)=0;
        end
        
    end
end

image_bitmap = bitmap;
end