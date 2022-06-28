function [triangle_matrix,triangle_polygon,xside] = generate_triangle(coords,xdim,ydim,dim)

if length(xdim)~=length(ydim)
    disp('Bounding X/Y regions must be the same size');
    return
end

%Specify bottom left corner coordinates
x0=coords(1); y0=coords(2);
xside = linspace(x0,x0+(xdim),dim);
yside = linspace(y0,y0+ydim,dim);
xc = x0 + xdim/2;
yc = y0 + ydim/3;

triangle_polygon = nsidedpoly(3,'Center',[xc,yc],'sidelength',xdim);
triangle_matrix = zeros([dim,dim]);
size(triangle_matrix)

for ii=1:dim
    for jj=1:dim
        x=xside(ii); 
        y=yside(jj);
        
        vertX = triangle_polygon.Vertices(:,1);
        vertY = triangle_polygon.Vertices(:,2);
        [in_poly,on_poly] = inpolygon(x,y,vertX,vertY);
        
        if in_poly ==1 || on_poly == 1
            triangle_matrix(ii,jj) = 1; %Flip bit-sign if in polygon
        end
    end
end


%Plotting of grid
%[X,Y] = meshgrid(xside,yside);
figure(),
plot(triangle_polygon,'facecolor','b'); grid on;
axis tight;
%plot(X,Y,'color','k'); hold on;
%plot(Y,X,'color','k'); hold on;
%xlim([(xside(1)-0.1*max(abs(xside))), xside(end)+0.1*max(abs(xside))]);
%xlim([(yside(1)-0.1*max(abs(yside))), yside(end)+0.1*max(abs(yside))]);

%% Testing bitmap update

% for kk1=1:dim
%     for kk2=1:dim
%     
%     xx = xside(kk1);
%     yy = yside(kk2);
%     [in_poly,on_poly] = inpolygon(xx,yy,vertX,vertY);
%     figure(123)
%     plot(X,Y,'color','k'); hold on;
%     plot(Y,X,'color','k'); hold on;
%     if in_poly ==1 || on_poly == 1
%         plot(xx,yy,'o','markersize',12,'linewidth',2,'color','b'); hold on; 
%     else
%         plot(xx,yy,'x','markersize',12,'linewidth',2,'color','r'); hold on; 
%     end
%     pause(0.005);
%     
%     end
% end
    


end