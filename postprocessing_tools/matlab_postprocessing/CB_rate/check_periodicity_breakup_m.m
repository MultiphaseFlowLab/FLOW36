function [ index ] = check_periodicity_breakup(xj,yj,zj,i,j,xn,xn1,corr )
% check all periodicity for breakup events

global mass
global step
global Lx
global Ly

% check original position and periodic positions
pos=[xj,yj,zj];
index1=check_closest(pos,xn1,j);
index=find(corr==index1);
% check volume correspondence
Vn=xn(index,7)*mass(i+step+1)/mass(i+1);
Vn1=xn1(j,7)+xn1(index1,7);
delta_V=abs((Vn-Vn1)/Vn);

% check all possible conditions at boundaries
% 1+3 possible positions at the corners, 1+1 at edges (original position included)
if(xj<=0.2*Lx && yj<=0.2*Ly)
    % lower left corner
    
    pos=[xj+Lx,yj,zj];
    index_1=check_closest(pos,xn1,j);
    index_01=find(corr==index_1);
    Vn=xn(index_01,7)*mass(i+step+1)/mass(i+1); %#ok<*FNDSB>
    Vn1=xn1(j,7)+xn1(index_1,7);
    if(abs((Vn-Vn1)/Vn)<delta_V)
        index=index_1;
        delta_V=abs((Vn-Vn1)/Vn);
    end
    
    pos=[xj,yj+Ly,zj];
    index_1=check_closest(pos,xn1,j);
    index_01=find(corr==index_1);
    Vn=xn(index_01,7)*mass(i+step+1)/mass(i+1);
    Vn1=xn1(j,7)+xn1(index_1,7);
    if(abs((Vn-Vn1)/Vn)<delta_V)
        index=index_1;
        delta_V=abs((Vn-Vn1)/Vn);
    end
    
    pos=[xj+Lx,yj+Ly,zj];
    index_1=check_closest(pos,xn1,j);
    index_01=find(corr==index_1);
    Vn=xn(index_01,7)*mass(i+step+1)/mass(i+1);
    Vn1=xn1(j,7)+xn1(index_1,7);
    if(abs((Vn-Vn1)/Vn)<delta_V)
        index=index_1;
%         delta_V=abs((Vn-Vn1)/Vn);
    end
    
%     maps(index,2)=j;
    
elseif(xj<=0.2*Lx && yj>=0.8*Ly)
    % lower right corner
    
    pos=[xj+Lx,yj,zj];
    index_1=check_closest(pos,xn1,j);
    index_01=find(corr==index_1);
    Vn=xn(index_01,7)*mass(i+step+1)/mass(i+1);
    Vn1=xn1(j,7)+xn1(index_1,7);
    if(abs((Vn-Vn1)/Vn)<delta_V)
        index=index_1;
        delta_V=abs((Vn-Vn1)/Vn);
    end
    
    pos=[xj,yj-Ly,zj];
    index_1=check_closest(pos,xn1,j);
    index_01=find(corr==index_1);
    Vn=xn(index_01,7)*mass(i+step+1)/mass(i+1);
    Vn1=xn1(j,7)+xn1(index_1,7);
    if(abs((Vn-Vn1)/Vn)<delta_V)
        index=index_1;
        delta_V=abs((Vn-Vn1)/Vn);
    end
    
    pos=[xj+Lx,yj-Ly,zj];
    index_1=check_closest(pos,xn1,j);
    index_01=find(corr==index_1);
    Vn=xn(index_01,7)*mass(i+step+1)/mass(i+1);
    Vn1=xn1(j,7)+xn1(index_1,7);
    if(abs((Vn-Vn1)/Vn)<delta_V)
        index=index_1;
%         delta_V=abs((Vn-Vn1)/Vn);
    end
    
%     maps(index,2)=j;
    
elseif(xj>=0.8*Lx && yj>=0.8*Ly)
    % upper right corner
    
    pos=[xj-Lx,yj,zj];
    index_1=check_closest(pos,xn1,j);
    index_01=find(corr==index_1);
    Vn=xn(index_01,7)*mass(i+step+1)/mass(i+1);
    Vn1=xn1(j,7)+xn1(index_1,7);
    if(abs((Vn-Vn1)/Vn)<delta_V)
        index=index_1;
        delta_V=abs((Vn-Vn1)/Vn);
    end
    
    pos=[xj,yj-Ly,zj];
    index_1=check_closest(pos,xn1,j);
    index_01=find(corr==index_1);
    Vn=xn(index_01,7)*mass(i+step+1)/mass(i+1);
    Vn1=xn1(j,7)+xn1(index_1,7);
    if(abs((Vn-Vn1)/Vn)<delta_V)
        index=index_1;
        delta_V=abs((Vn-Vn1)/Vn);
    end
    
    pos=[xj-Lx,yj-Ly,zj];
    index_1=check_closest(pos,xn1,j);
    index_01=find(corr==index_1);
    Vn=xn(index_01,7)*mass(i+step+1)/mass(i+1);
    Vn1=xn1(j,7)+xn1(index_1,7);
    if(abs((Vn-Vn1)/Vn)<delta_V)
        index=index_1;
%         delta_V=abs((Vn-Vn1)/Vn);
    end
    
%     maps(index,2)=j;
    
elseif(xj>=0.8*Lx && yj<=0.2*Ly)
    % upper left corner
    
    pos=[xj-Lx,yj,zj];
    index_1=check_closest(pos,xn1,j);
    index_01=find(corr==index_1);
    Vn=xn(index_01,7)*mass(i+step+1)/mass(i+1);
    Vn1=xn1(j,7)+xn1(index_1,7);
    if(abs((Vn-Vn1)/Vn)<delta_V)
        index=index_1;
        delta_V=abs((Vn-Vn1)/Vn);
    end
    
    pos=[xj,yj+Ly,zj];
    index_1=check_closest(pos,xn1,j);
    index_01=find(corr==index_1);
    Vn=xn(index_01,7)*mass(i+step+1)/mass(i+1);
    Vn1=xn1(j,7)+xn1(index_1,7);
    if(abs((Vn-Vn1)/Vn)<delta_V)
        index=index_1;
        delta_V=abs((Vn-Vn1)/Vn);
    end
    
    pos=[xj-Lx,yj+Ly,zj];
    index_1=check_closest(pos,xn1,j);
    index_01=find(corr==index_1);
    Vn=xn(index_01,7)*mass(i+step+1)/mass(i+1);
    Vn1=xn1(j,7)+xn1(index_1,7);
    if(abs((Vn-Vn1)/Vn)<delta_V)
        index=index_1;
%         delta_V=abs((Vn-Vn1)/Vn);
    end
    
%     maps(index,2)=j;
    
    % not in a corner, check for edges
elseif(xj<0.2*Lx)
    % lower edge
    
    pos=[xj+Lx,yj,zj];
    index_1=check_closest(pos,xn1,j);
    index_01=find(corr==index_1);
    Vn=xn(index_01,7)*mass(i+step+1)/mass(i+1);
    Vn1=xn1(j,7)+xn1(index_1,7);
    if(abs((Vn-Vn1)/Vn)<delta_V)
        index=index_1;
%         delta_V=abs((Vn-Vn1)/Vn);
    end
    
%     maps(index,2)=j;
    
elseif(yj<0.2*Ly)
    % left edge
    
    pos=[xj,yj+Ly,zj];
    index_1=check_closest(pos,xn1,j);
    index_01=find(corr==index_1);
    index_01=index_01(1);
    Vn=xn(index_01,7)*mass(i+step+1)/mass(i+1);
    Vn1=xn1(j,7)+xn1(index_1,7);
    if(abs((Vn-Vn1)/Vn)<delta_V)
        index=index_1;
%         delta_V=abs((Vn-Vn1)/Vn);
    end
    
%     maps(index,2)=j;
    
elseif(xj>0.8*Lx)
    % upper edge
    
    pos=[xj-Lx,yj,zj];
    index_1=check_closest(pos,xn1,j);
    index_01=find(corr==index_1);
    Vn=xn(index_01,7)*mass(i+step+1)/mass(i+1);
    Vn1=xn1(j,7)+xn1(index_1,7);
    if(abs((Vn-Vn1)/Vn)<delta_V)
        index=index_1;
%         delta_V=abs((Vn-Vn1)/Vn);
    end
    
%     maps(index,2)=j;
    
elseif(yj>0.8*Ly)
    % right edge
    
    pos=[xj,yj-Ly,zj];
    index_1=check_closest(pos,xn1,j);
    index_01=find(corr==index_1);
    Vn=xn(index_01,7)*mass(i+step+1)/mass(i+1);
    Vn1=xn1(j,7)+xn1(index_1,7);
    if(abs((Vn-Vn1)/Vn)<delta_V)
        index=index_1;
%         delta_V=abs((Vn-Vn1)/Vn);
    end
    
%     maps(index,2)=j;
    
end



end

