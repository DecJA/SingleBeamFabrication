function [sphere_3d_bitmap] = generate_3d_sphere_grid(diam,dim,minRad,coords)

xgrid=linspace(-diam/2,diam/2,dim);% + coords(1);
ygrid=xgrid;% + coords(2);
theta2 = linspace(pi/2,-pi/2,dim);

Rmax = diam/2;
Rheight = Rmax.*sin(theta2);

%%
sphere_3d_bitmap = zeros([dim dim length(Rheight)]);
coords = [0,0];
center = [coords(1),coords(2)];
%center=[(max(xgrid)+min(xgrid))/2+coords(1),(max(ygrid)+min(ygrid))/2+coords(2)];

for kk=1:length(Rheight)
    counter = 1;
    zPos = Rmax.*sin(theta2(kk));
    Rdisk = Rmax.*cos(theta2(kk));
    if Rdisk<minRad
        continue
    else
        for ii=1:dim
            for jj=1:dim

                posx=xgrid(ii);
                posy=ygrid(jj);
               % r = atan2(coords,[posx,posy]);
                r= sqrt((posx-center(1))^2 + ((posy-center(2))^2));

                if r<Rdisk && r>minRad 
                    sphere_3d_bitmap(ii,jj,kk) =1;
                    xloc(kk,counter) = xgrid(ii);
                    yloc(kk,counter) = ygrid(jj);
                    zloc(kk,counter) = zPos;
                    counter = counter+1;
                end
            end
        end
    end
end


noLayers = length(xloc(:,1));%length(find(abs(Rheight(:))>minRad))
for tt=1:noLayers
    A = find(xloc(tt,:)>0);
    B = find(yloc(tt,:)>0);
    if length(A)==0 || length(B)==0
        continue
    else
        figure(111),
        plot3(coords(1)+xloc(tt,:),coords(2)+yloc(tt,:),ones(size(xloc(tt,:))).*Rheight(tt),...
               '.','color','k'); hold on; grid on; grid minor;
    end
end
