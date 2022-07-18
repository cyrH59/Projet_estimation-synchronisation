function [dt,c] = preambuledetect(seqapp,rn,N)
% N representent le nombre de d pour lequel on veut tester la valeur 
c=zeros(1,N);

for k=1:N
 c(k)=sum(rn(k:(100+k-1)).*conj(seqapp));
end

[~,dt]=max(c);

end