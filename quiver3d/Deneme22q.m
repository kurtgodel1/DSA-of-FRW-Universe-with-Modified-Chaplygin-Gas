clear

syms x1 x2 x3 x4 x real

m=0.27;
% eqn1=-1-x3-3*x2+x1^2-x1*x3+x4==0;
% eqn2=x1*x3/m-x2*(2*x3-4-x1)==0;
% eqn3=x1*x3/m-2*x3*(x3-2)==0;
% eqn4=-2*x3*x4+x1*x4==0;


eqn1=-1-x3-3*x2+x1^2-x1*x3==0;
eqn2=x1*x3/m-x2*(2*x3-4-x1)==0;
eqn3=x1*x3/m-2*x3*(x3-2)==0;
eqn4=-2*x3*x4+x1*x4==0;


eqns=[eqn1 eqn2 eqn3];

[solx1,solx2,solx3,parameters,conditions]=solve(eqns,[x1 x2 x3],'ReturnConditions',true);
% 
% [solx1,solx2,solx3]=solve(eqns,[x1 x2 x3])
%clear


% for omega=-1:1:1
% for x1=-7:0.5:7
%     for x2=-7:0.5:7
%         u=-1-(omega-x1-x2)-3*x2+x1^2-x1*(omega-x1-x2);
%         v=x1*(omega-x1-x2)-x2*(2*(omega-x1-x2)-4-x1);
%         t=x1*(omega-x1-x2)-2*(omega-x1-x2)*((omega-x1-x2)-2);
%         
%         if u==0 & v==0 & t==0
%             plot3(x1,x2,(omega-x1-x2),'.','MarkerSize',20,'Color','k')
%         end
%         
%         l=sqrt(u^2+v^2+t^2);
%         
%         quiver3(x1,x2,(omega-x1-x2),u/l,v/l,t/l,.5,'MaxHeadSize',50,'Color','b')
%         hold on
%     end
% end
% end


[X Y Z]=meshgrid(-6:2:6,-6:2:6,0);
U=-1-Z-3*Y+X.^2-X.*Z;
V=X.*Z-Y.*(2*Z-4-X);
T=X.*Z-2*Z.*(Z-2);
L=sqrt(U.^2+V.^2+T.^2);

for ii=1:length(X)
quiver3D([X(:,ii),Y(:,ii),Z(:,ii)],[U(:,ii)./L(:,ii),V(:,ii)./L(:,ii),T(:,ii)./L(:,ii)],'b');
hold on
end
% 
% plot3([-1 1 0 -4],[0 0 -1 5],[0 0 2 0],'o')
% hold on

f=@(t,x)[-1-x(3)-3*x(2)+x(1)^2-x(1)*x(3);
        x(1)*x(3)/m-x(2)*(2*x(3)-4-x(1));
        x(1)*x(3)/m-2*x(3)*(x(3)-2)];
 
    
% solx1=double(solx1);
% solx2=double(solx2);
% solx3=double(solx3);
% 
% plot3(solx1,solx2,solx3,'.','MarkerSize',20);
% hold on
% 
% interval=[-0.1 0.1];
% differ=20;
%  for ii=1:length(solx1)
%      for jj=interval
%          for kk=interval
%              for ll=interval
%                   [t ya]=ode45(f,[0 20],[solx1(ii)-jj,solx2(ii)-kk,solx3(ii)-ll]);
%                   plot3(ya(:,1),ya(:,2),ya(:,3))
%                   hold on
% %                   len=floor(length(ya(:,1))/2);
% %                   xval1=ya(len,1);
% %                   xval2=ya(len+differ,1);
% %                   yval1=ya(len,2);
% %                   yval2=ya(len+differ,2);
% %                   zval1=ya(len,3);
% %                   zval2=ya(len+differ,3);
% %                   u=xval2-xval1;
% %                   v=yval2-yval1;
% %                   tt=zval2-zval1;
% %                   
% %                   quiver3(xval1,yval1,zval1,u,v,tt,'MaxHeadSize',1000);
% %                   [t ya]=ode45(f,[0 -20],[solx1(ii)-jj,solx2(ii)-kk,solx3(ii)-ll]);
% %                   plot3(ya(:,1),ya(:,2),ya(:,3))
%              end
%          end
%      end
%  end
% for ii=0.01:0.01:0.05
%     for jj=[.99 1.01 -.99 -1.01 -.01 .01 -4.01]
%     [t ya]=ode45(f,[0 20],[jj ii 1-ii-jj]);
%     plot3(ya(:,1),ya(:,2),ya(:,3))
%     hold on
%     [t ya]=ode45(f,[0 -20],[jj ii 1-ii-jj]);
%     plot3(ya(:,1),ya(:,2),ya(:,3))
%     end
% end
% 
% for ii=[-4.01 -3.99]
%     [t ya]=ode45(f,[0 20],[ii 1-ii 0]);
%     plot3(ya(:,1),ya(:,2),ya(:,3))
%     hold on
%     [t ya]=ode45(f,[0 -20],[ii 1-ii 0]);
%     plot3(ya(:,1),ya(:,2),ya(:,3))
% end    

% for ii=-6:6
%     for jj=-5:5
%         [t ya]=ode45(f,[0 20],[ii jj 1-ii-jj]);
%         plot3(ya(:,1),ya(:,2),ya(:,3))
%         hold on
%         [t ya]=ode45(f,[0 -20],[ii jj 1-ii-jj]);
%         plot3(ya(:,1),ya(:,2),ya(:,3))
%     end    
% end

% for ii=-6:2
%     for jj=-4:6
%         [t ya]=ode45(f,[0 20],[ii jj 0]);
%         plot3(ya(:,1),ya(:,2),ya(:,3))
%         hold on
%         [t ya]=ode45(f,[0 -20],[ii jj 0]);
%         plot3(ya(:,1),ya(:,2),ya(:,3))
%     end    
% end
xlabel('x1')
ylabel('x2')
zlabel('x3')

grid on
axis([-10 10 -10 10 -10 10])