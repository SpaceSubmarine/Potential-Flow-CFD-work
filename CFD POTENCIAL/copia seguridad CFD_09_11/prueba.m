clear
clc
close all;
tic
 
%% Input-DATA==============================================================
r=3.2;        %Radius of the cylinder (m)
v=30;        %Velocity (m/s)
dens=1025;  %density of sea water (kg/m3)
H=10;        %Height (m)
L=20;       %lenght (m)
tol=0.4;    %tolerance (hay problemas en funcion de la tolerancia que se escribe)
vOut=v;
maxIter=1e5;
maxdifer=1e-6;
 
%% Mesh====================================================================
dx=tol;
dy=tol;
M=round(H/dy);
N=round(L/dx);
X=linspace(0,L,N+2);
Y=linspace(H,0,M+2)';
[x,y]=meshgrid(X,Y);
 
%% Mfluidity===============================================================
Mfluid=ones(size(x));
Mdens=((ones(size(Mfluid)))*dens);
 
 
%% Central Point of Radius=================================================
%Mfluid(round((M+2)/2),round((N+2)/2))=0;
py=round((M+2)/2);
px=round((N+2)/2);
 
bP=(zeros(size(Mfluid)));
 
%% Cylinder================================================================
for j=(px-r/dx):(px+r/dx)
    for i=(py-r/dy):(py+r/dy)
        if (sqrt((x(i,j)-x(py,px)).^2+(y(i,j)-y(py,px)).^2))<r
            Mfluid(i,j)=0;
            bP(i,j)=v*H/2;
        end
    end
end
 
%% a_Coefficients==========================================================
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
 
 
Mdens=Mdens.*Mfluid;
 
%if(Mfluid == 1);
%% Esto arrglea los coeficientes en el interior del cilindro===============
 
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
 
 
 
%% GAUSS-SEIDEL METHOD=====================================================
 
PSI=(y*v).*Mfluid;
Mdens=Mdens.*Mfluid;
difer=inf;
Iter=1;
bP(:,end)=vOut*Y;
difer = Inf;
Iter = 0;
PSIOld=PSI;
tol = 0.0001;
 
while difer>maxdifer && Iter<maxIter
    
    PSIold=PSI;
    for i=2:M+1
        for j=2:N+1
            if(Mfluid(i,j)==1)
                ae(i,j)=((dPE/((dPe/(dens/Mdens(i,j)))+...
                    (dEe/(dens/Mdens(i,j+1)))))*(dy/dPE));
                
                ind=find(isinf(ae));
                ae(ind)=1;
                ind=find(isnan(ae));
                ae(ind)=1;
                
                aw(i,j)=(dPW/((dPw/(dens/Mdens(i,j)))+...
                    (dWw/(dens/Mdens(i,j-1)))))*(dy/dPW);
                
                ind=find(isinf(aw));
                aw(ind)=1;
                ind=find(isnan(aw));
                aw(ind)=1;
                
                as(i,j)=(dPS/((dPs/(dens/Mdens(i,j)))+...
                    (dSs/(dens/Mdens(i+1,j)))))*(dx/dPS);
                
                ind=find(isinf(as));
                as(ind)=1;
                ind=find(isnan(as));
                as(ind)=1;
                
                an(i,j)=(dPN/((dPn/(dens/Mdens(i,j)))+...
                    (dNn/(dens/Mdens(i-1,j)))))*(dx/dPN);
                
                ind=find(isinf(an));
                an(ind)=1;
                ind=find(isnan(an));
                an(ind)=1;
                
                ap(i,j)=ae(i,j)+aw(i,j)+as(i,j)+an(i,j);
                
                ind=find(isinf(ap));
                ap(ind)=1;
                ind=find(isnan(ap));
                ap(ind)=1;
                
            end
            
            PSI(i,j) = (ae(i,j)*PSI(i,j+1)+aw(i,j)*PSI(i,j-1)+...
                an(i,j)*PSI(i-1,j)+as(i,j)*PSI(i+1,j)+bP(i,j))/ap(i,j);
            
            
        end
    end
    Iter=Iter+1;
    difer=max(max(abs(PSIold-PSI)));
end
PSI(:,end)= PSI(:,end-1);
 
%% VELOCITY================================================================
vxP = v * Mfluid;
vyP = zeros(size(Mfluid));
vP = v * Mfluid;
vxn=ones(size(Mfluid));
vxs=ones(size(Mfluid));
vye=ones(size(Mfluid));
vyw=ones(size(Mfluid));
 
for i=2:M+1
    for j=2:N+1
        if(Mfluid(i,j)==1)
            
            vxn=an(i,j)*(PSI(i-1,j)-PSI(i,j))/dPN;
            vxs=as(i,j)*(PSI(i,j)-PSI(i+1,j))/dPS;
            vye=-ae(i,j)*(PSI(i,j+1)-PSI(i,j))/dPE;
            vyw=-aw(i,j)*(PSI(i,j)-PSI(i,j-1))/dPW;
            
            vxP(i,j)=(vxn+vxs)/2;
            vyP(i,j)=(vye+vyw)/2;
            
            vP(i,j)=sqrt(vxP(i,j)^2+vyP(i,j)^2);
            
            vP(i,j)=sqrt(vxP(i,j)^2+vyP(i,j)^2);
            
        end
    end
end
 
PSI=PSI.*Mfluid;
 
%% Plots===================================================================
 
% Stream function
% figure(1)
% contour(X,Y,PSI)
% title('Stream function')
% xlabel('X (m)')
% ylabel('Y (m)')
% xlim([0 L])
% ylim([0 H])
% title('Stream function')
% xlabel('X(m)')
% ylabel('Y(m)')
%
% figure(2)
% imagesc(X,Y,PSI)
% title('Stream function')
% xlabel('X (m)')
% ylabel('Y (m)')
% xlim([0 L])
% ylim([0 H])
% title('Stream function')
% xlabel('X(m)')
% ylabel('Y(m)')
%
% Velocity distribution
figure(3);
quiver(X,Y,vxP,vyP)
xlim([0 L])
ylim([0 H])
title('Velocity distribution')
xlabel('X')
ylabel('Y')
 
% % Mesh distribution
figure(4)
heatmap(Mdens);
title('Nodes distribution');
xlabel('X');
ylabel('Y');
%% Velocity_Plot===========================================================
% Create figure
figure1 = figure('WindowState','maximized',...
    'Colormap',[0 0 0;0.015873015873016 0.015873015873016 0.015873015873016;0.031746031746032 0.031746031746032 0.031746031746032;0.047619047619048 0.047619047619048 0.047619047619048;0.063492063492064 0.063492063492064 0.063492063492064;0.079365079365079 0.079365079365079 0.079365079365079;0.095238095238095 0.095238095238095 0.095238095238095;0.111111111111111 0.111111111111111 0.111111111111111;0.126984126984127 0.126984126984127 0.126984126984127;0.142857142857143 0.142857142857143 0.142857142857143;0.158730158730159 0.158730158730159 0.158730158730159;0.174603174603175 0.174603174603175 0.174603174603175;0.19047619047619 0.19047619047619 0.19047619047619;0.206349206349206 0.206349206349206 0.206349206349206;0.222222222222222 0.222222222222222 0.222222222222222;0.238095238095238 0.238095238095238 0.238095238095238;0.253968253968254 0.253968253968254 0.253968253968254;0.26984126984127 0.26984126984127 0.26984126984127;0.285714285714286 0.285714285714286 0.285714285714286;0.301587301587302 0.301587301587302 0.301587301587302;0.317460317460317 0.317460317460317 0.317460317460317;0.333333333333333 0.333333333333333 0.333333333333333;0.349206349206349 0.349206349206349 0.349206349206349;0.365079365079365 0.365079365079365 0.365079365079365;0.380952380952381 0.380952380952381 0.380952380952381;0.396825396825397 0.396825396825397 0.396825396825397;0.412698412698413 0.412698412698413 0.412698412698413;0.428571428571429 0.428571428571429 0.428571428571429;0.444444444444444 0.444444444444444 0.444444444444444;0.46031746031746 0.46031746031746 0.46031746031746;0.476190476190476 0.476190476190476 0.476190476190476;0.492063492063492 0.492063492063492 0.492063492063492;0.507936507936508 0.507936507936508 0.507936507936508;0.523809523809524 0.523809523809524 0.523809523809524;0.53968253968254 0.53968253968254 0.53968253968254;0.555555555555556 0.555555555555556 0.555555555555556;0.571428571428571 0.571428571428571 0.571428571428571;0.587301587301587 0.587301587301587 0.587301587301587;0.603174603174603 0.603174603174603 0.603174603174603;0.619047619047619 0.619047619047619 0.619047619047619;0.634920634920635 0.634920634920635 0.634920634920635;0.650793650793651 0.650793650793651 0.650793650793651;0.666666666666667 0.666666666666667 0.666666666666667;0.682539682539683 0.682539682539683 0.682539682539683;0.698412698412698 0.698412698412698 0.698412698412698;0.714285714285714 0.714285714285714 0.714285714285714;0.73015873015873 0.73015873015873 0.73015873015873;0.746031746031746 0.746031746031746 0.746031746031746;0.761904761904762 0.761904761904762 0.761904761904762;0.777777777777778 0.777777777777778 0.777777777777778;0.793650793650794 0.793650793650794 0.793650793650794;0.80952380952381 0.80952380952381 0.80952380952381;0.825396825396825 0.825396825396825 0.825396825396825;0.841269841269841 0.841269841269841 0.841269841269841;0.857142857142857 0.857142857142857 0.857142857142857;0.873015873015873 0.873015873015873 0.873015873015873;0.888888888888889 0.888888888888889 0.888888888888889;0.904761904761905 0.904761904761905 0.904761904761905;0.920634920634921 0.920634920634921 0.920634920634921;0.936507936507937 0.936507936507937 0.936507936507937;0.952380952380952 0.952380952380952 0.952380952380952;0.968253968253968 0.968253968253968 0.968253968253968;0.984126984126984 0.984126984126984 0.984126984126984;1 1 1]);
 
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
 
% Create image
image(vP,'Parent',axes1,'CDataMapping','scaled');
 
% Create zlabel
zlabel('ZLabel');
 
% Create ylabel
ylabel('Y Position');
 
% Create xlabel
xlabel({'X Position',''});
 
% Create title
title('Velocity Field');
 
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[1.41469223184944 48.6874195045767]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[-6.51912873111979 33.8586895341725]);
box(axes1,'on');
axis(axes1,'ij');
% Set the remaining axes properties
set(axes1,'CLim',[10.1334427406547 101.334427406547],'Colormap',...
    [0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 0.9375;0.125 1 0.875;0.1875 1 0.8125;0.25 1 0.75;0.3125 1 0.6875;0.375 1 0.625;0.4375 1 0.5625;0.5 1 0.5;0.5625 1 0.4375;0.625 1 0.375;0.6875 1 0.3125;0.75 1 0.25;0.8125 1 0.1875;0.875 1 0.125;0.9375 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0;0.5 0 0],...
    'DataAspectRatio',[1 1 1],'Layer','top','PlotBoxAspectRatio',...
    [23.6363636363636 20.1889091326461 1]);
% Create colorbar
colorbar(axes1);
axis equal;
hold off;
 
%% Streamlines_Plot========================================================
% Create figure
figure1 = figure('WindowState','maximized');
 
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
 
% Create contour
contour(PSI,'LineWidth',0.01,'LineColor','none',...
    'LevelStep',1,...
    'Fill','on');
 
% Create ylabel
ylabel('Y(m)');
 
% Create xlabel
xlabel('X(m)');
 
% Create title
title('Stream function');
 
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 20]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0 10]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','Color',[0 0 0],'Colormap',...
    [0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 0.9375;0.125 1 0.875;0.1875 1 0.8125;0.25 1 0.75;0.3125 1 0.6875;0.375 1 0.625;0.4375 1 0.5625;0.5 1 0.5;0.5625 1 0.4375;0.625 1 0.375;0.6875 1 0.3125;0.75 1 0.25;0.8125 1 0.1875;0.875 1 0.125;0.9375 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0;0.5 0 0],...
    'Layer','top');
axis equal;
hold off;
 
toc

