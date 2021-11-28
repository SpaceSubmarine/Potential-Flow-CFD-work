%% My velocity


vxP = v * Mfluid;
vyP = zeros(size(Mfluid));
vP = v * Mfluid;

for i=2:M+1
    for j=2:N+1
        if(Mfluid(i,j)==1)
            
            vxn(i,j)=(dPN/((dPn/(dens/Mdens(i,j)))+...
                (dNn/(dens/Mdens(i-1,j)))))...
                *((PSI(i-1,j)-PSI(i,j))/dPN);
            
            vxs(i,j)=(dPS/((dPs/(dens/Mdens(i,j)))+...
                (dSs/(dens/Mdens(i+1,j)))))...
                *((PSI(i,j)-PSI(i-1,j))/dPS);
            
            vye(i,j)=-(dPE/((dPe/(dens/Mdens(i,j)))+...
                (dEe/(dens/Mdens(i,j+1)))))...
                *((PSI(i,j+1)-PSI(i,j))/dPE);
            
            vyw(i,j)=-(dPW/((dPw/(dens/Mdens(i,j)))+...
                (dWw/(dens/Mdens(i,j-1)))))...
                *((PSI(i,j)-PSI(i,j-1))/dPW);
            
            vxP(i,j)=(vxn(i,j)+vxs(i,j))/2;
            vyP(i,j)=(vye(i,j)+vyw(i,j))/2;
            
            vP(i,j)=sqrt(vxP(i,j)^2+vyP(i,j)^2);
            
        end
    end
end
