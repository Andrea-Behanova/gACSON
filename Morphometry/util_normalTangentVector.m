function nTangVec = util_normalTangentVector(curve,disp) 

x = curve(:,1);
y = curve(:,2);

dx = gradient(x);
dy = gradient(y);
ds = sqrt(dx.^2 + dy.^2);
nTangVec = [dx./ds, dy./ds];

if disp == 1 && size(curve,2)==2
    plot(y,x,'Color',[0,0.66,0],'LineWidth',2)
    hold on
    quiver(y,x,nTangVec(:,2),nTangVec(:,1),0.5,'b');
    grid on;
    xlabel('x'); ylabel('y');
end

if size(curve,2)==3
    z = curve(:,3);
    dz = gradient(z);
    ds = sqrt(dx.^2 + dy.^2 + dz.^2);
    nTangVec = [dx./ds, dy./ds, dz./ds];
    
    if disp == 1
        plot3(y,x,z,'Color',[0,0.66,0],'LineWidth',2)
        hold on
        %quiver3(y,x,z,nTangVec(:,2),nTangVec(:,1),nTangVec(:,3),0.5,'b');
        %grid on;
        %xlabel('x'); ylabel('y'); zlabel('z')
        axis equal
    end
end




