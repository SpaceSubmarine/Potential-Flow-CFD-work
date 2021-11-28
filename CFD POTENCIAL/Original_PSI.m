%% PSI=====================================================================
while difer>maxdifer && Iter<maxIter
    
    PSIold=PSI;
    for i=2:M+1
        for j=2:N+1
            
            PSI(i,j) = (ae(i,j)*PSI(i,j+1)+aw(i,j)*PSI(i,j-1)+...
                an(i,j)*PSI(i-1,j)+as(i,j)*PSI(i+1,j)+bP(i,j))/ap(i,j);
            
        end
    end
    Iter=Iter+1;
    difer=max(max(abs(PSIold-PSI)));
end
PSI(:,end)= PSI(:,end-1);
PSI=PSI.*Mfluid;
%PSI=flipud(PSI);