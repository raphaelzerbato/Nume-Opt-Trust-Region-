clear ''ALL'' 
close all 
clc
%% figure start
figure(1) 
title('path to the minimum')
hold on
%% line search

x= [-50;-100];
delta=1;
image_1=f_2(x);
[delta,dkvect] = trust_method_NM1(x,delta);
%eps = norm(dkvect);

%%%%%graph
figure (1)
supx =abs(x(1,1))+1;
supy =abs(x(2,1))+1;
step = 1;
xabs = zeros(1,round(supx*step*2));
yord=zeros(1,round(supy*step*2));

xabs(1,1) = -(supx);
yord(1,1)= -(supy);

for i=1:2*supx
 xabs(1,i+1) = xabs(1,i)+step;
end

for i=1:2*supy
    yord(1,i+1)=yord(1,i)+step;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION TO CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y]=meshgrid(xabs,yord);
F = X.^4+Y.^2+9.*(X.*Y).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


levels= 0:(image_1/5):image_1;
[C,h] =contour(X,Y,F,round(levels));%'LineColor','blue','ShowText','on');
clabel(C,h,'LabelSpacing',1000,'FontSize',15);
plot(x(1,1),x(2,1),'or','MarkerSize',5)
drawnow
%%% FIN GRAPH%%%%%%%
%%% new point 

newcord = x+(delta*dkvect);
image_2 = f_2(newcord);
[delta,dkvect,H] =trust_method_NM1(x,delta);
eps = norm(dkvect);
plot(newcord(1,1),newcord(2,1),'or','MarkerSize',5)
%% algorithme
while eps>0.01
    
    x = newcord;
    plot(x(1,1),x(2,1),'or','MarkerSize',5)
    [delta,dkvect,H] = trust_method_NM1(x,delta);
    eps = norm(dkvect);
   
    %%% new point
    newcord =x+(delta*dkvect);
    image_2 = f_2(newcord);
end
hold off
