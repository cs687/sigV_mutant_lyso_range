a=[];
for i=1:length(loc)+1; 
    if i==1; 
        a=len(1:loc(i),1);
    elseif i==length(loc)+1; 
        a=[a;a(end)+[len(loc(i-1)+1:end,1)]];
    else; 
        a=[a;a(end)+[len(loc(i-1)+1:loc(i),1)]]; 
    end
end