function visualize(u,v,p,space,time,res,FV)
% This function plots the results for velocits and pressure
close all

figure(1)
dp = mean2(p(:,end/2))-mean(p(:,end));

usample = u(:,end-2);
ysample = space.Y(:,end-1);
uanalytical = -15*(ysample.*(ysample-space.h));
plot(usample,ysample,'LineWidth',2)
hold on
plot(uanalytical,ysample,'LineWidth',2)
legend("Numerical solution","Analytical solution")
xlabel("u Velocity [m/s]")
ylabel("Y [m]")
set(gca, 'FontSize',20)
% 
% x = linspace(0,space.l,space.dimX);
% for j=2:space.dimX
%     vector = u(:,j);
%     k = find(vector>0.1);
%     delta_exp(j) = space.Y(k(end),1);
% end
% 
% figure(2)
% pcolor(space.X,space.Y,u);
% shading interp
% c = colorbar;
% c.Label.String = 'U Velocity [m/s]';
% xlabel('X [m]')
% ylabel('Y [m]')
% hold on
% plot(space.X(2,1:end),delta_exp,'Color','k')
% set(gca, 'FontSize',20)
%hold on
%plot(x,5*sqrt((FV.nu*x)/(0.15)))

% % Create residual plot
% figure(1)
% plot(res,'LineWidth',3)
% set(gca, 'YScale', 'log')
% set(gca, 'FontSize',20)
% xlabel('Number of Iterations','LineWidth',3)
% ylabel('Residual','LineWidth',3)
% 
% % Create u velocity plot in 3D and 2D
% figure(2)
% %subplot(2,1,1)
% %surf(space.X,space.Y,u)
% %xlabel('X [m]')
% %ylabel('Y [m]')
% %zlabel('Velocity u [m/2]')
% %subplot(2,1,2)
% pcolor(space.X,space.Y,u);
% shading interp
% c = colorbar;
% c.Label.String = 'U Velocity [m/s]';
% xlabel('X [m]')
% ylabel('Y [m]')
% set(gca, 'FontSize',20)
% 
% % Create v velocity plot
% figure(3)
% %subplot(2,1,1)
% %surf(space.X,space.Y,v)
% %xlabel('X [m]')
% %ylabel('Y [m]')
% %zlabel('Velocity v [m/s]')
% %subplot(2,1,2)
% pcolor(space.X,space.Y,v)
% shading interp
% c = colorbar;
% c.Label.String = 'V Velocity [m/s]';
% xlabel('X [m]')
% ylabel('Y [m]')
% set(gca, 'FontSize',20)  
% 
% % Create pressure plot
% %subplot(2,1,1)
% %surf(space.X,space.Y,p)
% %xlabel('X [m]')
% %ylabel('Y [m]')
% %zlabel('Pressure [Pa]')
% %subplot(2,1,2)
% pcolor(space.X,space.Y,p)
% shading interp
% c = colorbar;
% c.Label.String = 'Pressure [Pa]';
% xlabel('X [m]')
% ylabel('Y [m]')
% set(gca, 'FontSize',20)
% 
% % Reduce matrix size for vector field visualization
% unew=u;
% vnew=v;
% for i=1:space.dimY
%     for j=1:space.dimX
%         if(mod(i,4)~=0 || mod(j,4)~=0)
%             unew(i,j) = 0;
%             vnew(i,j) = 0;
%         end
%     end
% end
% 
% % Create vector plot
% 
% figure(5)
% quiver(space.X,space.Y,unew,vnew,5)
% xlim([0 space.l])
% ylim([0 space.h])
% hold on
% pos = 0.4;
% r = 0.1;
% formfunction = @(xnorm) real(sqrt(r^2-(xnorm-pos)^2));
% fplot(formfunction,'Color','k')
% set(gca, 'FontSize',20)
%  
% % Create streamline plot
% figure(6)
% startx = 1.0*ones(30,1);
% starty = linspace(0,space.h,30);
% streamline(space.X,space.Y,u,v,startx,starty,[0.001,1e6])
% hold on
% fplot(formfunction,'Color','k')
% xlim([0,space.l])
% ylim([0,space.h])
% set(gca, 'FontSize',20)

end