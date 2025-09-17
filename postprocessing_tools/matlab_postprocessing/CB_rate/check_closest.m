function [ index ] = check_closest(xn,xn1,j)
% find closest position to xn in xn1

delta=sqrt((xn(1)-xn1(:,1)).^2+(xn(2)-xn1(:,2)).^2+(xn(3)-xn1(:,3)).^2);
% exclude droplet itself
delta(j)=NaN;

% find closest droplet
[~,index]=min(delta);

end

