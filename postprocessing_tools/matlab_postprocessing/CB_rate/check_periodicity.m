function [ indmin ] = check_periodicity(xj,yj,zj,xn1)
% check periodicity (initial check, no volume check, smallest distance)

global Lx
global Ly

% check original position and periodic positions
pos=[xj,yj,zj];
r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
[~,indmin]=min(r);
rmin=r(indmin);

% check all possible conditions at boundaries
% 1+3 possible positions at the corners, 1+1 at edges (original position included)
if(xj<=0.2*Lx && yj<=0.2*Ly)
    % lower left corner
    
    pos=[xj+Lx,yj,zj];
    r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
    [~,indmin_1]=min(r);
    rmin_1=r(indmin_1);
    if(rmin_1<rmin)
        rmin=rmin_1;
        indmin=indmin_1;
    end
    
    pos=[xj,yj+Ly,zj];
    r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
    [~,indmin_1]=min(r);
    rmin_1=r(indmin_1);
    if(rmin_1<rmin)
        rmin=rmin_1;
        indmin=indmin_1;
    end
    
    pos=[xj+Lx,yj+Ly,zj];
    r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
    [~,indmin_1]=min(r);
    rmin_1=r(indmin_1);
    if(rmin_1<rmin)
        %         rmin=rmin_1;
        indmin=indmin_1;
    end
    
elseif(xj<=0.2*Lx && yj>=0.8*Ly)
    % lower right corner
    
    pos=[xj+Lx,yj,zj];
    r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
    [~,indmin_1]=min(r);
    rmin_1=r(indmin_1);
    if(rmin_1<rmin)
        rmin=rmin_1;
        indmin=indmin_1;
    end
    
    pos=[xj,yj-Ly,zj];
    r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
    [~,indmin_1]=min(r);
    rmin_1=r(indmin_1);
    if(rmin_1<rmin)
        rmin=rmin_1;
        indmin=indmin_1;
    end
    
    pos=[xj+Lx,yj-Ly,zj];
    r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
    [~,indmin_1]=min(r);
    rmin_1=r(indmin_1);
    if(rmin_1<rmin)
        %         rmin=rmin_1;
        indmin=indmin_1;
    end
    
elseif(xj>=0.8*Lx && yj>=0.8*Ly)
    % upper right corner
    
    pos=[xj-Lx,yj,zj];
    r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
    [~,indmin_1]=min(r);
    rmin_1=r(indmin_1);
    if(rmin_1<rmin)
        rmin=rmin_1;
        indmin=indmin_1;
    end
    
    pos=[xj,yj-Ly,zj];
    r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
    [~,indmin_1]=min(r);
    rmin_1=r(indmin_1);
    if(rmin_1<rmin)
        rmin=rmin_1;
        indmin=indmin_1;
    end
    
    pos=[xj-Lx,yj-Ly,zj];
    r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
    [~,indmin_1]=min(r);
    rmin_1=r(indmin_1);
    if(rmin_1<rmin)
        %         rmin=rmin_1;
        indmin=indmin_1;
    end
    
elseif(xj>=0.8*Lx && yj<=0.2*Ly)
    % upper left corner
    
    pos=[xj-Lx,yj,zj];
    r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
    [~,indmin_1]=min(r);
    rmin_1=r(indmin_1);
    if(rmin_1<rmin)
        rmin=rmin_1;
        indmin=indmin_1;
    end
    
    pos=[xj,yj+Ly,zj];
    r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
    [~,indmin_1]=min(r);
    rmin_1=r(indmin_1);
    if(rmin_1<rmin)
        rmin=rmin_1;
        indmin=indmin_1;
    end
    
    pos=[xj-Lx,yj+Ly,zj];
    r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
    [~,indmin_1]=min(r);
    rmin_1=r(indmin_1);
    if(rmin_1<rmin)
        %         rmin=rmin_1;
        indmin=indmin_1;
    end
    
    % not in a corner, check for edges
elseif(xj<0.2*Lx)
    % lower edge
    pos=[xj+Lx,yj,zj];
    r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
    [~,indmin_1]=min(r);
    rmin_1=r(indmin_1);
    if(rmin_1<rmin)
        %         rmin=rmin_1;
        indmin=indmin_1;
    end
    
elseif(yj<0.2*Ly)
    % left edge
    pos=[xj,yj+Ly,zj];
    r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
    [~,indmin_1]=min(r);
    rmin_1=r(indmin_1);
    if(rmin_1<rmin)
        %         rmin=rmin_1;
        indmin=indmin_1;
    end
    
elseif(xj>0.8*Lx)
    % upper edge
    pos=[xj-Lx,yj,zj];
    r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
    [~,indmin_1]=min(r);
    rmin_1=r(indmin_1);
    if(rmin_1<rmin)
        %         rmin=rmin_1;
        indmin=indmin_1;
    end
    
elseif(yj>0.8*Ly)
    % right edge
    pos=[xj,yj-Ly,zj];
    r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
    [~,indmin_1]=min(r);
    rmin_1=r(indmin_1);
    if(rmin_1<rmin)
        %         rmin=rmin_1;
        indmin=indmin_1;
    end
    
    
end

end

