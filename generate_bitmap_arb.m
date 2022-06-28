function [image_bitmap] = generate_bitmap_arb(coords,xlength,ylength,dim,image)
addpath('R:\OMG\declan\PrintingSummaryProgram\Runfiles\PrintPatterns');
x0 = coords(1);
y0 = coords(2);

%Scale ratio of points
max_side = find(max([xlength,ylength]));
if max_side == 1
    xdim = dim;
    ydim = (xlength/ylength)*dim;
elseif max_side == 2
    ydim = dim;
    xdim = (ylength/xlength)*dim;
end

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