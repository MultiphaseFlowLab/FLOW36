function [ Af ] = filter_sum( A,w )
% filter data

l=int16(length(A(:,1))/w)-1;

Af=zeros(l,2);

j=0;
for i=1:w:length(A(:,1))-w+1
    j=j+1;
    Af(j,2)=sum(A(i:i+w-1,2));
    Af(j,1)=mean(A(i:i+w-1,1));
end


end

