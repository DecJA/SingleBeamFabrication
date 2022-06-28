function [voxels,bitmap] = generate_sphere_even_1hole(diam,spacing,hdiam)
% diam = 4e-6;
% hdiam = 1e-6;
% spacing = 1064e-9/10;
radius = diam/2;
hradius = hdiam/2;
numr = ceil(radius/spacing) + 0.5; %even dipole range 
dim = numr - 0.5;
rrange = (-numr:numr)*spacing;

% Generate the voxel grid
[xx,yy,zz] = meshgrid(rrange, rrange, rrange);

%calculate hole centeres - hard code to start
[r, t, p] = ott.utils.xyz2rtp(xx, yy, zz); %Full sphere coordinates
c1 = -radius/2; 
c2 = radius/2;

b = r<radius;

%Remove first hole
xx2 = xx(b); yy2 = yy(b); zz2 = zz(b); %regular sphere
[~,rlayer] = cart2pol(xx2+c1,yy2);%translate sphere center over hole
bRemove = (rlayer>hradius);

voxels = [xx2(bRemove), yy2(bRemove), zz2(bRemove)].';

figure(),
plot3(xx2(bRemove),yy2(bRemove),zz2(bRemove),'.','color','k');

%Generate bitmap for voxels
bitmap = zeros([dim dim dim]);

noVox = size(voxels(1,:),2);

for n=1:noVox
    r= voxels(:,n)';
    
    diffX = abs(rrange-r(1));
    diffY = abs(rrange-r(2));
    diffZ = abs(rrange-r(3));
    
    [~,idx] = min(diffX);
    [~,idy] = min(diffY);
    [~,idz] = min(diffZ);
    
    bitmap(idx,idy,idz) = 1;

end





