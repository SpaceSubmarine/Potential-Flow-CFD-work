for i=2:M+1
    
    for j=2:N+1
        
        ae(i,j)=((dPE/((dPe/(dens/Mdens(i,j)))+...
            (dEe/(dens/Mdens(i,j+1)))))*(dy/dPE));
        %{
        ind=find(isinf(ae));
        ae(ind)=0;
        ind=find(isnan(ae));
        ae(ind)=0;
            %}
            aw(i,j)=(dPW/((dPw/(dens/Mdens(i,j)))+...
                (dWw/(dens/Mdens(i,j-1)))))*(dy/dPW);
            %{
        ind=find(isinf(aw));
        aw(ind)=0;
        ind=find(isnan(aw));
        aw(ind)=0;
                %}
                as(i,j)=(dPS/((dPs/(dens/Mdens(i,j)))+...
                    (dSs/(dens/Mdens(i+1,j)))))*(dx/dPS);
                %{
        ind=find(isinf(as));
        as(ind)=0;
        ind=find(isnan(as));
        as(ind)=0;
                    %}
                    an(i,j)=(dPN/((dPn/(dens/Mdens(i,j)))+...
                        (dNn/(dens/Mdens(i-1,j)))))*(dx/dPN);
                    %{
        ind=find(isinf(an));
        an(ind)=0;
        ind=find(isnan(an));
        an(ind)=0;
                        %}
                        ap(i,j)=ae(i,j)+aw(i,j)+as(i,j)+an(i,j);
                        %{
        ind=find(isinf(ap));
        ap(ind)=1;
        ind=find(isnan(ap));
        ap(ind)=1;
                        %}
    end
    
end
