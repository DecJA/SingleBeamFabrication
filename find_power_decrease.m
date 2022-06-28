function [powerSurface,Xq,Yq] = find_power_decrease(xl,yl)

%Experimental polynomial fit
funcFit = @(p00,p10,p01,p20,p11,p02,p30,p21,p12,p03,x,y) (p00 + p10.*x + p01.*y + p20.*x.^2 + p11.*x.*y + p02.*y.^2 + p30.*x.^3 + p21.*x.^2.*y + p12.*x.*y.^2 + p03.*y.^3);
p00 =      0.9975;
p10 =        -193;
p01 =      -107.6;
p20 =  -2.834e+06;
p11 =  -1.279e+05;
p02 =  -3.211e+06;
p30 =   7.627e+11;
p21 =   2.835e+11;
p12 =   2.649e+11;
p03 =   6.936e+09;


dxl = 100e-9;
xq = (min(xl):dxl:max(xl));
yq = (min(yl):dxl:max(yl));
[X,Y] = meshgrid(xl,yl);
surfRough = funcFit(p00,p10,p01,p20,p11,p02,p30,p21,p12,p03,X,Y);
[Xq,Yq] = meshgrid(xq,yq);
powerSurface = interp2(X,Y,surfRough,Xq,Yq);

end