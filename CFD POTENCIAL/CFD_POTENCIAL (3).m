clear
clc
tic

%Input-DATA------------------------------------------------------------
r=3.2;        %Radius of the cylinder (m)
v=8;        %Velocity (m/s)
dens=1025;  %density of sea water (kg/m3)
H=10;        %Height (m)
L=20;       %lenght (m)
tol=0.4;    %tolerance (hay problemas en funcion de la tolerancia que se escribe)
vOut=v;
maxIter=1e5;
maxdifer=1e-6;

%Mesh------------------------------------------------------------------
dx=tol;
dy=tol;
M=round(H/dy);
N=round(L/dx);
X=linspace(0,L,N+2);
Y=linspace(H,0,M+2)';
[x,y]=meshgrid(X,Y);

%Mfluidity---------------------------------------------------------------
Mfluid=ones(size(x));
Mdens=((ones(size(Mfluid)))*dens);


%------------Central Point of Radius-----------------------------------
%Mfluid(round((M+2)/2),round((N+2)/2))=0;
py=round((M+2)/2);
px=round((N+2)/2);

bP=(zeros(size(Mfluid)));

%=========================Cylinder===================================
for j=(px-r/dx):(px+r/dx)
   for i=(py-r/dy):(py+r/dy)
      if (sqrt((x(i,j)-x(py,px)).^2+(y(i,j)-y(py,px)).^2))<r
         Mfluid(i,j)=0;
         bP(i,j)=v*y(i,j)/2;
      end
   end
end

%=============================a_Coefficients==============================
ae=(ones(size(Mfluid)));
aw=(ones(size(Mfluid)));
as=(ones(size(Mfluid)));            
an=(ones(size(Mfluid)));            
ap=(ones(size(Mfluid))); 
     

%conditions
bP(:,1)=v*Y;
an(:,1)=0;
as(:,1)=0;
ae(:,1)=0;
aw(:,1)=0;
an(:,end)=0;
as(:,end)=0;
ae(:,end)=0;
bP(1,:)=H*v;

dPE=dx;
dPW=dx;
dPS=dy;
dPN=dy;
dPe=dx/2;
dPw=dx/2;
dPs=dy/2;
dPn=dy/2;
dEe=dx/2;
dWw=dx/2;
dSs=dy/2;
dNn=dy/2;

PSI=(y*v).*Mfluid;
Mdens=Mdens.*Mfluid;
difer=inf;
Iter=1;
bP(:,end)=vOut*Y;

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
    
%Esto arrglea los coeficientes en el interior del cilindro
for j=(px-r/dx):(px+r/dx)
   for i=(py-r/dy):(py+r/dy)
      if (sqrt((x(i,j)-x(py,px)).^2+(y(i,j)-y(py,px)).^2))<r
         ap(i,j)=1;
         ae(i,j)=0;
         aw(i,j)=0;
         as(i,j)=0;
         an(i,j)=0;
      end
   end
end
    
%=============================PSI==============================    
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

PSI=PSI.*Mfluid;



%ALGO SALE MAL
%============================VELOCITY==================================
vxn=(ones(size(Mfluid)));
vye=(ones(size(Mfluid)));
vxs=(ones(size(Mfluid)));
vyw=(ones(size(Mfluid)));
vxp=(ones(size(Mfluid)));
vyp=(ones(size(Mfluid)));
vP=(ones(size(Mfluid)));
for i=2:M+1
    for j=2:N+1
     
        vxn(i,j)=(dPN/((dPn/(dens/Mdens(i,j)))+...
            (dNn/(dens/Mdens(i-1,j)))))...
            *(PSI(i-1,j)-PSI(i,j))/dPN;
        
        vxs(i,j)=(dPS/((dPs/(dens/Mdens(i,j)))+...
            (dSs/(dens/Mdens(i+1,j)))))...
            *(PSI(i+1,j)-PSI(i,j))/dPS;
        
        vye(i,j)=(dPE/((dPe/(dens/Mdens(i,j)))+...
            (dEe/(dens/Mdens(i,j+1)))))...
            *(PSI(i,j+1)-PSI(i,j))/dPE;
        
        vyw(i,j)=(dPW/((dPw/(dens/Mdens(i,j)))+...
            (dWw/(dens/Mdens(i,j-1)))))...
            *(PSI(i,j-1)-PSI(i,j))/dPW;
        
        vxp=(vxn-vxs)/2;
        vyp=(vye-vyw)/2;
        
        vP(i,j)=sqrt((vxp(i,j)^2)+(vyp(i,j)^2));
        
        if (sqrt((x(i,j)-x(py,px)).^2+(y(i,j)-y(py,px)).^2))<r
         vxn(i,j)=0;
         vxs(i,j)=0;
         vye(i,j)=0;
         vyw(i,j)=0;
         vP(i,j)=0;
        end
        
    end
end

    
%-------------Plots----------------------------------------------------
%figure;heatmap(PSI);
%figure;heatmap(Mfluid);
figure;imagesc(PSI);
%figure;heatmap(ap);
%figure;heatmap(Mdens);
%hold on;
%figure;contour(vy);
%figure;quiver(vx,vy);
%hold off
%figure;heatmap(vy);
%figure;heatmap(vx);
%-------------Plots----------------------------------------------------

toc