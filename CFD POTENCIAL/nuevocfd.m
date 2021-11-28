clear;
close all; 
clc;
tic

%% INPUTS==================================================================

v=60;                       %(m/s)  
H=10;                        %(m)
L=20;                        %(m)
T = 288.15;                 %kelvin
P = 101325;                 %Pa
R = 286.9;                  %[J/K/mol]
dens = P/R/T;                %[Kg/m^3]
gamma = 1.41;
cP = gamma*R/(gamma-1);     %J/kg/K

r = L/10;                   %(m)
               %(m)

maxIter=1e4;
tol=0.4;    %tolerance (hay problemas en funcion de la tolerancia que se escribe)
%% Meshgrid
%{
N = 60;
M = 50;

dx = L/N;
dy = H/M;

% The boundary nodes are not at the same length of the other ones
% so they must be add after the mesh creation

xp = linspace(0, L, N+2);
yp = linspace(0, H, M+2);
[X,Y] = meshgrid(xp,yp);

matrix = ones(size(X));
%}

% Mesh====================================================================
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

%% Cylinder
%Central Point of Radius=================================================
%Mfluid(round((M+2)/2),round((N+2)/2))=0;
py=round((M+2)/2);
px=round((N+2)/2);

PSI= (y*v).*(Mfluid);
for j=(px-r/dx):(px+r/dx)
    for i=(py-r/dy):(py+r/dy)
        if (sqrt((x(i,j)-x(py,px)).^2+(y(i,j)-y(py,px)).^2))<r
            Mfluid(i,j)=0;
            PSI(i,j)=v*H/2;
        end
    end
end


%% Setting initial parameters

Ti = T * Mfluid;
Pi = P * Mfluid; 
Di = dens * Mfluid;

%% Setting distances

dPE = dx;
dPW = dx;
dPS = dy;
dPN = dy;
dPe = dx/2;
dPw = dx/2;
dPs = dy/2;
dPn = dy/2;
dEe = dx/2;
dWw = dx/2;
dSs = dy/2;
dNn = dy/2;

%% Setting coefficients

ae=(ones(size(Mfluid)));
aw=(ones(size(Mfluid)));
as=(ones(size(Mfluid)));            
an=(ones(size(Mfluid)));            
ap=(ones(size(Mfluid))); 

%% Boundary conditions

an(:,1) = 0;
as(:,1) = 0;
ae(:,1) = 0;
aw(:,1) = 0;
an(:,end) = 0;
as(:,end) = 0;
ae(:,end) = 0;

%% Gauss-Seidel

global_error = Inf;
n_iter = 0;
PSIOld = PSI;
tol = 0.0001;

while ((global_error>tol) && (n_iter<maxIter))
    error = zeros(M+2,N+2);
    for i=2:M+1
        for j=2:N+1
            
            PSIOld = PSI(i,j);
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
            
            error(i,j) = abs(PSIOld-PSI(i,j));
            if error (i,j) > tol
                error(i,j) = 1;
            end
        end
    end
    
    global_error =sum(sum(error));
    n_iter = n_iter+1;
end

PSI(:,end)= PSI(:,end-1);
PSI= PSI.*Mfluid;

%% Calculate velocities

vxP = v * Mfluid;
vyP = zeros(M+2,N+2);
vP = v * Mfluid;

for i=2:M+1
        for j=2:N+1                  
            if(Mfluid(i,j)==1)
                
                vxn(i,j)=an(i,j)*(PSI(i-1,j)-PSI(i,j))/dPN;
                vxs(i,j)=as(i,j)*(PSI(i+1,j)-PSI(i-1,j))/dPS;
                vye(i,j)=-ae(i,j)*(PSI(i,j+1)-PSI(i,j))/dPE;
                vyw(i,j)=-aw(i,j)*(PSI(i,j-1)-PSI(i,j))/dPW;
    
                vxP(i,j)=(vxn(i,j)+vxs(i,j))/2;
                vyP(i,j)=(vye(i,j)+vyw(i,j))/2;
  
                vP(i,j)=sqrt(vxP(i,j)^2+vyP(i,j)^2);
                
            end
        end
end

%% Calculate final parameters

V_0 = v * ones(M+2,N+2);
v_0_square = V_0.*V_0;
vP_square = vP.*vP; 

Tp = (Ti + (v_0_square - vP_square)/(2*cP));
Pp = Pi.*(Tp./Ti).^(gamma/(gamma-1));
dens_P = Pp./(R*Tp);

%% Plots

% Stream function

    figure(1)
    contour(X,Y,PSI)
    title('Stream function')
    xlabel('X (m)')
    ylabel('Y (m)')
    xlim([0 L])
    ylim([0 H])
    title('Stream function')
    xlabel('X(m)')
    ylabel('Y(m)')

% Velocity distribution

    figure(2);
    quiver(X,Y,vxP,vyP)
    xlim([0 L])
    ylim([0 H])
    title('Velocity distribution')
    xlabel('X(m)')
    ylabel('Y(m)')
    
% Pressure distribution

    figure(3)
    pcolor(X,Y,Pp); 
    shading interp
    axis([0, L, 0, H])
    daspect([1, 1, 1])
    colorbar
    title('Pressure distribution (Pa)')
    xlabel('X(m)')
    ylabel('Y(m)')
    
% Temperature distribution

    figure(4);
    hold on
    pcolor(X,Y,Tp); 
    shading interp
    axis([0, L, 0, H])
    daspect([1, 1, 1])
    colorbar
    title('Temperature distribution (K)')
    xlabel('X(m)')
    ylabel('Y(m)')
hold off
% Density distribution

    figure(5)
    hold on
    pcolor(X,Y,dens_P); 
    shading interp
    axis([0, L, 0, H])
    daspect([1, 1, 1])
    colorbar
    title('Density distribution (Kg/m^3)')
    xlabel('X(m)')
    ylabel('Y(m)')
    hold off
% Mesh distribution

    figure(6) 
    plot(x,y,'.','Color','r')
    hold on
    title('Nodes distribution')
    xlabel('X(m)')
    ylabel('Y(m)')
    hold off
toc
