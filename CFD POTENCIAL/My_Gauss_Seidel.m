difer = Inf;
Iter = 0;
PSI = PSI_i;
tol = 0.0001;

while ((difer>tol) && (Iter<maxIter))
    error = zeros(M+2,N+2);
    for i=2:M+1
        for j=2:N+1
            
            PSIold = PSI(i,j);
            value = false;
            
            if (Mfluid(i,j)==1)
                
                ae(i,j)=(dPE/((dPe/(dens/Mdens(i,j)))+...
                    (dEe/(dens/Mdens(i,j+1)))))*(dy/dPE);
                
                if(ae(i,j)==Inf || isnan(ae(i,j)))
                    ae(i,j)=1;
                    value=true;
                end
                
                aw(i,j)=(dPW/((dPw/(dens/Mdens(i,j)))+...
                    (dWw/(dens/Mdens(i,j-1)))))*(dy/dPW);
                
                if(aw(i,j)==Inf || isnan(aw(i,j)))
                    aw(i,j)=1;
                    value=true;
                end
                
                as(i,j)=(dPS/((dPs/(dens/Mdens(i,j)))+...
                    (dSs/(dens/Mdens(i+1,j)))))*(dx/dPS);
                
                if(as(i,j)==Inf || isnan(as(i,j)))
                    as(i,j)=1;
                    value=true;
                end
                
                an(i,j)=(dPN/((dPn/(dens/Mdens(i,j)))+...
                    (dNn/(dens/Mdens(i-1,j)))))*(dx/dPN);
                
                if(an(i,j)==Inf || isnan(an(i,j)))
                    an(i,j)=1;
                    value=true;
                end
                
                if(value)
                    ap(i,j)=1;
                else
                    ap(i,j)=ae(i,j)+aw(i,j)+as(i,j)+an(i,j);
                end
                
                PSI(i,j)=(ae(i,j)*PSI(i,j+1)+aw(i,j)*PSI(i,j-1)+...
                    an(i,j)*PSI(i-1,j)+as(i,j)*PSI(i+1,j))/ap(i,j);
            end
            
            error(i,j) = abs(PSIold-PSI(i,j));
            if error (i,j) > tol
                error(i,j) = 1;
            end
        end
    end
    
    difer =sum(sum(error));
    Iter = Iter+1;
end

PSI(:,end)= PSI(:,end-1);
PSI_e = PSI.*Mfluid;