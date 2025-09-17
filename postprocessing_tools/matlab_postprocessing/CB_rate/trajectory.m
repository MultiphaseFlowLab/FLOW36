% PROBLEM: as it is now, considers only binary breakups

% try to match all droplets without volume check, quality is then reduced
% if no volume matching is found
% verifies better population balance, worse in volume matching

clear all
close all
clc

figure(1)
hold all
% 
% figure(2)
% hold all

global mass
global step
global Lx
global Ly

init=560000;
step=2500;
last=750000;

timespan=int16((last-init)/step);
B=zeros(timespan,2);
C=zeros(timespan,2);

Re=300;
Lx=4*pi*Re;
Ly=2*pi*Re;

% tolerance on volume loss [0;1]
voltol=0.2;

path='../../coalescence_breakup/output/';

% used to estimate mass losses
A=importdata(['../../results/','time_check.dat'],' ',1);
mass=A.data(:,4);


A=importdata([path,'drop_count.dat'],' ',1);
indstart=find(A.data(:,1)==init);
number=A.data(:,3);

Deltat=1.5e-5;
dt=step*Deltat*Re;

fprintf(' Iteration      N +     B -     C =   Bal |   N+1   Delta   Quality \n');

for i=init:step:last-step
    A=importdata([path,'center_timestep_',num2str(i,'%08i'),'.dat'],' ',1);
    xn=A.data(:,2:8);
    n=length(xn(:,1));
    

    xnest(1:n,1)=xn(:,1)+xn(:,4)*dt;
    xnest(1:n,2)=xn(:,2)+xn(:,5)*dt;
    xnest(1:n,3)=xn(:,3)+xn(:,6)*dt;

    % check for periodicity
    for ii=1:n
        while xnest(ii,1)>Lx && xnest(ii,1)>0
            xnest(ii,1)=xnest(ii,1)-Lx;
        end
        while xnest(ii,2)>Ly && xnest(ii,2)>0
            xnest(ii,2)=xnest(ii,2)-Ly;
        end
        while xnest(ii,1)<0
            xnest(ii,1)=xnest(ii,1)+Lx;
        end
        while xnest(ii,2)<0
            xnest(ii,2)=xnest(ii,2)+Ly;
        end
    end
    
    
    corr=zeros(n,1);
    
    A=importdata([path,'center_timestep_',num2str(i+step,'%08i'),'.dat'],' ',1);
    xn1=A.data(:,2:8);
    n1=length(xn1(:,1));
    
    for j=1:n
        xj=xnest(j,1);
        yj=xnest(j,2);
        zj=xnest(j,3);
        if((xj>0.2*Lx) && (xj<0.8*Lx) && (yj>0.2*Ly) && (yj<0.8*Ly))
            pos=[xj,yj,zj];
            r=sqrt((pos(1)-xn1(:,1)).^2+(pos(2)-xn1(:,2)).^2+(pos(3)-xn1(:,3)).^2);
            [~,indmin]=min(r);
            % assign indmin to current drop
            corr(j)=indmin;
        else
            indmin=check_periodicity(xj,yj,zj,xn1);
            corr(j)=indmin;
        end
        
       
        
    end
    
    % initialize quality index to n, every time no volume correspondence is
    % found, reduce quality by 1
    quality=n;
    
    % initialize breakup and coalescence counters
    breakup=0;
    coalescence=0;
    
    % mapping array, array indexing correspond to step n, index contained
    % in array correspond to step n+1
    maps=corr;
    maps(:,2)=0;
    
    % start from breakups, find indexes of xn1 missing in corr
    % only bunary breakups
    for j=1:n1
        if(sum(ismember(j,corr))==0)
            % breakup occurred, drop j has no original droplet
            % fprintf('breakup\n');
            breakup=breakup+1;
            
            xj=xn1(j,1);
            yj=xn1(j,2);
            zj=xn1(j,3);
            % check if close to boundaries
            if((xj>0.2*Lx) && (xj<0.8*Lx) && (yj>0.2*Ly) && (yj<0.8*Ly))
                % if in inner domain no periodicity issues expected
                pos=[xj,yj,zj];
                index1=check_closest(pos,xn1,j);
                index=find(corr==index1);
                maps(index,2)=j;
            else
                index=check_periodicity_breakup(xj,yj,zj,i,j,xn,xn1,corr);
                maps(index,2)=j;
            end
        end
    end
    
    
    % check for coalescences, only binary coalescences
    while isempty(corr)==0
        index=find(corr==corr(1));
        if(length(index)>1)
            coalescence=coalescence+1;
        end
        corr(index)=[];        
    end
    
    
    
    % check volume using maps array and eventually reduce quality
    while (isempty(maps)==0 && isempty(xn)==0)
        index=find(maps==maps(1));
        
        if(length(index)>1)
            % coalescence
            Vn=(xn(index(1),7)+xn(index(2),7))*mass(i+step+1)/mass(i+1);
            Vn1=xn1(maps(index(1),1),7);
            delta_V=abs((Vn-Vn1)/Vn);
            if(delta_V<voltol)
                % correspondence found, check on volume passed
            else
                quality=quality-1;
            end
        elseif(maps(index,2)==0)
            % traslation
            Vn=xn(index(1),7)*mass(i+step+1)/mass(i+1);
            Vn1=xn1(maps(index(1),1),7);
            delta_V=abs((Vn-Vn1)/Vn);
            if(delta_V<voltol)
                % correspondence found, check on volume passed
            else
                quality=quality-1;
            end
        else
            % breakup
            Vn=xn(index(1),7)*mass(i+step+1)/mass(i+1);
            Vn1=xn1(maps(index(1),1),7)+xn1(maps(index(1),2),7);
            delta_V=abs((Vn-Vn1)/Vn);
            if(delta_V<voltol)
                % correspondence found, check on volume passed
            else
                quality=quality-1;
            end
        end
        % clear already processed entries
        maps(index,:)=[];   
        xn(index,:)=[];
    end
    
    figure(1)
    plot(i,quality/n*100,'ko')
%     
%     figure(2)
%     plot(i,n,'ko')
%     plot(i,breakup,'ro')
%     plot(i,coalescence,'bo')
    
    B((i-init)/step+1,:)=[i,breakup];
    C((i-init)/step+1,:)=[i,coalescence];

    % print balance and quality
    fprintf('%10i  %5i + %5i - %5i = %5i | %5i   %5i     %5.1f\n',i,n,breakup,coalescence,n+breakup-coalescence,n1,n1-(n+breakup-coalescence),quality/n*100);
    
end

figure
hold all
plot(B(:,1)*Re*Deltat,B(:,2)./number(indstart:end-1),'ro')
plot(C(:,1)*Re*Deltat,C(:,2)./number(indstart:end-1),'bo')
% plot(B(:,1)*Re*Deltat,B(:,2),'ro')
% plot(C(:,1)*Re*Deltat,C(:,2),'bo')

% w=10;
% Bf=filter_sum(B,w);
% Cf=filter_sum(C,w);
% numf=filter_sum([B(:,1),number(1:end-1,1)],w);
% numf=numf(:,2);

% plot(Bf(:,1)*Re*Deltat,Bf(:,2)./numf,'r')
% plot(Cf(:,1)*Re*Deltat,Cf(:,2)./numf,'b')
% plot(Bf(:,1)*Re*Deltat,Bf(:,2),'r')
% plot(Cf(:,1)*Re*Deltat,Cf(:,2),'b')


% tBC=[Bf(:,1),Bf(:,2)./numf,Cf(:,2)./numf];
% tBC=[B(:,1),B(:,2)./number(1:end-1),C(:,2)./number(1:end-1)];
% tBC=[B(:,1),B(:,2)./number(1),C(:,2)./number(1)];


tBC=[B(:,1)*Re*Deltat,B(:,2),B(:,2)./number(indstart:end-1),C(:,2),C(:,2)./number(indstart:end-1)];

save('tBC','tBC');



