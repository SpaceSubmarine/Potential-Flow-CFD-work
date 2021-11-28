clear
clc
close all;
tic

%% Input-DATA==============================================================
r=3.2;      %Radius of the cylinder (m)
v=30;       %Velocity (m/s)
dens=1025;  %density of sea water (kg/m3)
H=10;       %Height (m)
L=20;       %lenght (m)
tol=0.4;    %tolerance 
vOut=v;     %Velocity (m/s)
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
            %bP(i,j)=200;
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
%% This corrects the coefficients inside the cylinder======================

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
            
        end
    end
end

PSI=PSI.*Mfluid;

%% Velocity distribution Plots=============================================
figure(3);
quiver(X,Y,vxP,vyP)
xlim([0 L])
ylim([0 H])
title('Velocity distribution')
xlabel('X')
ylabel('Y')

%% Mesh distribution=======================================================
figure(4)
heatmap(Mdens);
title('Nodes Distribution and Density');
xlabel('X');
ylabel('Y');

%% Velocity_Plot===========================================================
% Create figure
figure1 = figure('WindowState','maximized',...
    'Colormap',[0.041666666666667 0 0;0.083333333333333 0 0;0.125 0 0;0.166666666666667 0 0;0.208333333333333 0 0;0.25 0 0;0.291666666666667 0 0;0.333333333333333 0 0;0.375 0 0;0.416666666666667 0 0;0.458333333333333 0 0;0.5 0 0;0.541666666666667 0 0;0.583333333333333 0 0;0.625 0 0;0.666666666666667 0 0;0.708333333333333 0 0;0.75 0 0;0.791666666666667 0 0;0.833333333333333 0 0;0.875 0 0;0.916666666666667 0 0;0.958333333333333 0 0;1 0 0;1 0.041666666666667 0;1 0.083333333333333 0;1 0.125 0;1 0.166666666666667 0;1 0.208333333333333 0;1 0.25 0;1 0.291666666666667 0;1 0.333333333333333 0;1 0.375 0;1 0.416666666666667 0;1 0.458333333333333 0;1 0.5 0;1 0.541666666666667 0;1 0.583333333333333 0;1 0.625 0;1 0.666666666666667 0;1 0.708333333333333 0;1 0.75 0;1 0.791666666666667 0;1 0.833333333333333 0;1 0.875 0;1 0.916666666666667 0;1 0.958333333333333 0;1 1 0;1 1 0.0625;1 1 0.125;1 1 0.1875;1 1 0.25;1 1 0.3125;1 1 0.375;1 1 0.4375;1 1 0.5;1 1 0.5625;1 1 0.625;1 1 0.6875;1 1 0.75;1 1 0.8125;1 1 0.875;1 1 0.9375;1 1 1]);

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
title('Velocity Field (m/s)');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0.5 102.5]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0.499999999999998 52.5]);
box(axes1,'on');
axis(axes1,'ij');
% Set the remaining axes properties
set(axes1,'CLim',[10.1334427406547 101.334427406547],'Colormap',...
    [0.041666666666667 0 0;0.083333333333333 0 0;0.125 0 0;0.166666666666667 0 0;0.208333333333333 0 0;0.25 0 0;0.291666666666667 0 0;0.333333333333333 0 0;0.375 0 0;0.416666666666667 0 0;0.458333333333333 0 0;0.5 0 0;0.541666666666667 0 0;0.583333333333333 0 0;0.625 0 0;0.666666666666667 0 0;0.708333333333333 0 0;0.75 0 0;0.791666666666667 0 0;0.833333333333333 0 0;0.875 0 0;0.916666666666667 0 0;0.958333333333333 0 0;1 0 0;1 0.041666666666667 0;1 0.083333333333333 0;1 0.125 0;1 0.166666666666667 0;1 0.208333333333333 0;1 0.25 0;1 0.291666666666667 0;1 0.333333333333333 0;1 0.375 0;1 0.416666666666667 0;1 0.458333333333333 0;1 0.5 0;1 0.541666666666667 0;1 0.583333333333333 0;1 0.625 0;1 0.666666666666667 0;1 0.708333333333333 0;1 0.75 0;1 0.791666666666667 0;1 0.833333333333333 0;1 0.875 0;1 0.916666666666667 0;1 0.958333333333333 0;1 1 0;1 1 0.0625;1 1 0.125;1 1 0.1875;1 1 0.25;1 1 0.3125;1 1 0.375;1 1 0.4375;1 1 0.5;1 1 0.5625;1 1 0.625;1 1 0.6875;1 1 0.75;1 1 0.8125;1 1 0.875;1 1 0.9375;1 1 1],...
    'DataAspectRatio',[1 1 1],'GridColor','none','Layer','top',...
    'PlotBoxAspectRatio',[51 26 1]);
colorbar(axes1);
hold off;

%% Streamlines_Plot========================================================
%CREATEFIGURE1(PSI)
%  PSI:  contour z

%  Auto-generated by MATLAB on 19-Nov-2020 12:20:54

% Create figure
figure1 = figure('WindowState','maximized');

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create contour
contour(PSI,'LineWidth',0.01,'LineColor','none','LevelStep',1,...
    'Fill','on');

% Create ylabel
ylabel('Y (M) Nodes');

% Create xlabel
xlabel('X (N) Nodes');

% Create title
title('Stream Function');

box(axes1,'on');
axis(axes1,'tight');
axis(axes1,'ij');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','Color',[0 0 0],'Colormap',...
    [0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 0.9375;0.125 1 0.875;0.1875 1 0.8125;0.25 1 0.75;0.3125 1 0.6875;0.375 1 0.625;0.4375 1 0.5625;0.5 1 0.5;0.5625 1 0.4375;0.625 1 0.375;0.6875 1 0.3125;0.75 1 0.25;0.8125 1 0.1875;0.875 1 0.125;0.9375 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0;0.5 0 0],...
    'DataAspectRatio',[1 1 1],'Layer','top','PlotBoxAspectRatio',[50.5 25.5 1]);
% Create colorbar
colorbar(axes1);
axis equal;
hold off;

%% Streamlines_Plot2=======================================================
% Create figure
figure1 = figure('WindowState','maximized');

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create contour
contour(PSI,'LineWidth',0.01,'LevelStep',5);

% Create ylabel
ylabel('Y (M) Nodes');

% Create xlabel
xlabel('X (N) Nodes');

% Create title
title('Stream Function');

box(axes1,'on');
axis(axes1,'tight');
axis(axes1,'ij');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','Color',[0 0 0],'Colormap',...
    [0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 0.9375;0.125 1 0.875;0.1875 1 0.8125;0.25 1 0.75;0.3125 1 0.6875;0.375 1 0.625;0.4375 1 0.5625;0.5 1 0.5;0.5625 1 0.4375;0.625 1 0.375;0.6875 1 0.3125;0.75 1 0.25;0.8125 1 0.1875;0.875 1 0.125;0.9375 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0;0.5 0 0],...
    'DataAspectRatio',[1 1 1],'Layer','top','PlotBoxAspectRatio',[50.5 25.5 1]);
% Create colorbar
colorbar(axes1);
axis equal;
hold off;
%%
toc
