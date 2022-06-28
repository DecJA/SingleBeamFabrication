function [sphere_3d_bitmap,pointLoc,counter2] = generate_3d_sphere_grid_2_holes(diam,dim,rad1,seperation,coords,Rmin)
% diam = 5-6;
% dim =80;
% rad1 = 600e-9;
%seteration500-800e-9
 rad2 = rad1;

xgrid=linspace(-diam/2,diam/2,dim);% + coords(1);
ygrid=xgrid;% + coords(2);
theta2 = linspace(pi/2,-pi/2,dim);

Rmax = diam/2;
Rheight = Rmax.*sin(theta2);

%%
sphere_3d_bitmap = zeros([dim dim length(Rheight)]);
mainCenter = [0,0];

xloc = zeros([dim,length(Rheight)]);
yloc = xloc;

%Determine location for center of 2 small spots
if rad1>= diam/2
    warning('Overlappling sphere holes. Reduce radius!');
    return
end
extraRadius = diam - (2*rad1) - 0.1*(diam/2); %Allow 10% of sphere edge - CHECK!
holeLoc = extraRadius/2; %Equally spaced 

% r1 = [rad1 + holeLoc/4,0];
mustBePositive(diam/2-(rad1+seperation/2));
r1 = [rad1+seperation/2,0];
r2 = -r1;
maxDist = abs(coords) + rad1;
if abs(maxDist(1))>abs(Rmax) || abs(maxDist(2))>abs(Rmax)
    disp('Holes exceede fill circule diamater, change hole size/position');
end

xvec = linspace(-diam/2,diam/2,dim)+coords(1);
yvec = linspace(-diam/2,diam/2,dim) +coords(2);
zvec = linspace(-diam/2,diam/2,length(Rheight));

counter2 = 1;
for kk=1:length(Rheight)
    counter = 1;
    zPos = Rmax.*sin(theta2(kk));
    Rdisk = Rmax.*cos(theta2(kk));
    if Rdisk<Rmin
       continue
    else
        for ii=1:dim
            for jj=1:dim
                posx = xgrid(ii);
                posy=ygrid(jj);

                rd= sqrt((posx-r1(1))^2 + ((posy-r1(2))^2));
                rd2= sqrt((posx-r2(1))^2 + ((posy-r2(2))^2));
                r= sqrt((posx)^2 + ((posy)^2));
                if (r<Rdisk)
                    if rd>rad1 && rd2>rad2
                        sphere_3d_bitmap(ii,jj,kk) =1;
                        xloc(kk,counter) = xgrid(ii);
                        yloc(kk,counter) = ygrid(jj);
                        zloc(kk,counter) = zPos;
                        counter = counter+1;
                        counter2 = counter2+1;
                        pointLoc(:,counter2) = [xvec(ii),yvec(jj),zvec(kk)];
%                         pointLoc(1,ii) = xvec(ii);
%                         pointLoc(2,jj) = yvec(jj);
%                         pointLoc(3,kk) = zvec(kk);
                    end
                end
            end
        end
    end
end


noLayers = length(xloc(:,1));

slice = floor(length(Rheight)/2);
figure(),
imagesc(sphere_3d_bitmap(:,:,slice).'); colormap gray
title('Slice from middle of shape');

% for tt=1:noLayers
%     figure(1),
%     plot3(xloc(tt,:),yloc(tt,:),ones(size(xloc(tt,:))).*Rheight(tt),...
%            '.','color','k'); hold on; grid on; grid minor;
% end

end