function dene(g)

    for ii=0:0.1:1
        for jj=0:0.1:(1-ii)
            
        [X, Y, Z]=meshgrid(ii,jj,1-ii-jj);
        


        U = ((3*g-2).*X-2*Y+2-3*g).*X;
        V = 2*(1+1/2*(3*g-2).*X-Y).*Y;
        T = ((3*g-2)*X-2*Y).*Z;

        L = sqrt(U.^2+V.^2+T.^2);

        quiver3D([X(:),Y(:),Z(:)],[U(:),V(:),T(:)]./100,'r',0.94)

        hold on
        
        end
    end
    
    [X, Y, Z] = meshgrid(0:0.1:1);
    
    U = ((3*g-2).*X-2*Y+2-3*g).*X;
    V = 2*(1+1/2*(3*g-2).*X-Y).*Y;
    T = ((3*g-2)*X-2*Y).*Z;
    
    startx=[];
    starty=[];
    startz=[];
    
    if g<=0.1
        d=0.05;
        e=3*d;
    else
        d=0.01;
        e=d;
    end
    
    
    for ii=0:0.1:1
        for jj=0:d:e
            startx=[startx ii];
            starty=[starty jj];
            startz=[startz 1-ii-jj];
        end
    end
    
    h=streamline(X,Y,Z,U,V,T,startx,starty,startz);
    set(h,'Color','g');
    set(h,'LineWidth',1.5);
    
    syms x y z

    [a,b,c] = solve(((3*g-2)*x-2*y+2-3*g)*x, 2*(1+1/2*(3*g-2)*x-y)*y, ((3*g-2)*x-2*y)*z);

    plot3(a,b,c,'o')
    
    axis([0 1 0 1 0 1])


    xlabel('\Omega')
    ylabel('\Omega_{\Lambda}')
    zlabel('K')
    view(0,90)
    hold off
end