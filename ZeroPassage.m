function [min] = ZeroPassage(Ftau)
N=length(Ftau);
mins=[];
for i=1:N-1
    if sign(Ftau(i))~=sign(Ftau(i+1))
        mins=[mins,i+1];
    end
end
if length(mins)>1
    for i =1 :length(mins)
        if(Ftau(mins(i))-Ftau(mins(i)-1)>0) %% max dérivée seconde 
            min=mins(i);
        end
    end
else
min=mins;
end
end

