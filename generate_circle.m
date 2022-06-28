function [circle_array] = generate_circle(coords,diam,dim)

%Specify bottom left grid location start, 
%diam is circle circumference.

x0=coords(1);
y0=coords(2);
xside=linspace(x0,x0+diam,dim);
yside=linspace(y0,y0+diam,dim);
[X,Y] = meshgrid(xside,yside);

center=[(max(xside)+min(xside))/2,(max(yside)+min(yside))/2];
circle_array = zeros([dim,dim]);

for ii=1:dim
    for jj=1:dim
        x=xside(ii);
        y=yside(jj);
        
        loc = sqrt(abs(x-center(1))^2 + abs((y-center(2)))^2);
        
        if loc <= diam/2
            circle_array(ii,jj) = 1;
        end
%         
%         figure(123),
%         plot(X,Y,'color','k'); hold on;
%         plot(Y,X,'color','k'); hold on;
%         if circle_array(ii,jj) == 0
%             plot(x,y,'o','markersize',12,'linewidth',2,'color','b'); hold on;
%         else
%             plot(x,y,'x','markersize',12,'linewidth',2,'color','r'); hold on;
%         end
    end
end
   
end