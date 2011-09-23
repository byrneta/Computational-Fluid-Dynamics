%%  Convection-Diffusion PDE

%% Definitions

clc; clear all; close all;      % Clear Workspace
L = 1;				% Length of Channel, L (m)
T = 0.2;			% Time, T (s)
a = 2.5;			% Velocity, a (m/s)
alpha = 0.005;			% Acceleration, alpha (m^2/s)
Uo = 1.0;                       % BC u(0,t), Uo (m/s)

%% Analytical Solution
xa = 0:0.001:L;                 % Create Analytic X Column Vector, xa
xo = 0.2;			% Initial Position, xo (m)

Ua = 1-0.5*(1+erf((xa-xo-a*T)/(2*sqrt(alpha*T))))'; 

%% Initialization 

dx = 0.015;                     % Step Size in X, dx
dt = (0.05*dx^2)/alpha;         % Step Size in T (Calculated), dt
d = alpha*dt/dx^2;              % Fourier Coefficient, d
c = a*dt/dx;                    % Courant Coefficient, c
Rec = a*dx/alpha;               % Cell Reynolds Number, Rec
x = [0:dx:L]';                  % Create X Column Vector, x
t = [0:dt:T];                   % Create Time Step Row Vector, t
numT = length(t);               % Number of Time Steps, numT
numX = length(x);               % Number of X-Direction Mesh Points, numX


U = zeros(numX,numT);           % Create Velocity Matrix

U(1,:) = Uo;                    % Apply BC: u(0,t) = Uo;
U(numX,:) = 0;                  % Apply BC: u(L,t) = 0;

for i=1:numX                    % Apply IC's
    if x(i)<0.2
       U(i,1) = 1.0;            % If x<0.2 ==> u(x,0)=1.0
    end
    if x(i)==0.2
        U(i,1) = 0.5;           % If x=0.2 ==> u(x,0)=0.5
    end
    if x(i)>0.2
        U(i,1) = 0.0;           % If x>0.2 ==> u(x,0)=0.0
    end
end

%% Finite Volume - FTCS Convection / FTCS Diffusion

Uvftcs=U;
for j=2:numT
    for i=2:numX-1
        Uvftcs(i,j) = Uvftcs(i,j-1) - (c/2)*( Uvftcs(i+1,j-1) - Uvftcs(i-1,j-1) ) + d*( Uvftcs(i+1,j-1) - 2*Uvftcs(i,j-1) + Uvftcs(i-1,j-1) );
    end
end

%% Finite Volume - First Order Upwind Convection / FTCS Diffusion

Uvuf=U;
for j=2:numT
    for i=2:numX-1
        Uvuf(i,j) = Uvuf(i,j-1) - c*( Uvuf(i,j-1) - Uvuf(i-1,j-1) ) + d*( Uvuf(i+1,j-1) - 2*Uvuf(i,j-1) + Uvuf(i-1,j-1) );
    end
end



%% Finite Difference - FTCS Convection / FTCS Diffusion

Uftcs=U;
for j=2:numT
    for i=2:numX-1
        Uftcs(i,j) = Uftcs(i,j-1) - (c/2)*( Uftcs(i+1,j-1) - Uftcs(i-1,j-1) ) + d*( Uftcs(i+1,j-1) - 2*Uftcs(i,j-1) + Uftcs(i-1,j-1) );
    end
end

%% Finite Difference - First Order Upwind Convection / FTCS Diffusion

Uuf=U;
for j=2:numT
    for i=2:numX-1
        Uuf(i,j) = Uuf(i,j-1) - c*( Uuf(i,j-1) - Uuf(i-1,j-1) ) + d*( Uuf(i+1,j-1) - 2*Uuf(i,j-1) + Uuf(i-1,j-1) );
    end
end

%% Finite Difference - Lax-Wendroff Convection / FTCS Diffusion

Ulwf=U;
for j=2:numT
    for i=2:numX-1
        Ulwf(i,j) = Ulwf(i,j-1) - (c/2)*( Ulwf(i+1,j-1) - Ulwf(i-1,j-1) ) + (c^2/2)*( Ulwf(i+1,j-1) - 2*Ulwf(i,j-1) + Ulwf(i-1,j-1) ) + d*( Ulwf(i+1,j-1) - 2*Ulwf(i,j-1) + Ulwf(i-1,j-1) );
    end
end

%% Finite Difference - MacCormack Convection / FTCS Diffusion

Umf=U;
Umfs=Umf;
%for j=2:numT
%    for i=2:numX-1
%        Umfs(i,j) = Umf(i,j-1) - c*( Umf(i+1,j-1) - Umf(i,j-1) ) + d*( Umf(i+1,j-1) - 2*Umf(i,j-1) + Umf(i-1,j-1) );
%        Umf(i,j) = 0.5*( Umf(i,j-1) + Umfs(i,j-1) ) - c*( Umfs(i,j) - Umfs(i-1,j) ) + d*( Umfs(i+1,j-1) - 2*Umfs(i,j-1) + Umfs(i-1,j-1) );
%    end
%end
for j=1:numT-1
    for i=2:numX-1
        Umfs(i,j) = Umf(i,j) - c*( Umf(i+1,j) - Umf(i,j) ) + d*( Umf(i+1,j) - 2*Umf(i,j) + Umf(i-1,j) );
    end
    for i=2:numX-1
        Umf(i,j+1) = 0.5*( Umf(i,j) + Umfs(i,j)  - c*( Umfs(i,j) - Umfs(i-1,j) )) + d*( Umfs(i+1,j) - 2*Umfs(i,j) + Umfs(i-1,j) );
    end
end


%% Plot Results
figure;
hold on
plot(x,Uftcs(:,numT),'ro--');
plot(x,Uuf(:,numT),'m+--');
plot(x,Ulwf(:,numT),'g*--');
plot(x,Umf(:,numT),'kx--');
plot(xa,Ua,'-');
ylabel('U (m/s)');
xlabel('X (m)');
legend('FTCS/FTCS ','Upwind/FTCS','Lax-Wendroff/FTCS','MacCormack/FTCS','Exact');
title(['Analytical vs Finite Difference Schemes at T=',num2str(T),'s , \Deltax=',num2str(dx),'m , Rec=',num2str(Rec)])

figure;
hold on
plot(x,Uvftcs(:,numT),'ro--');
plot(x,Uvuf(:,numT),'m+--');
plot(xa,Ua,'-');
ylabel('U (m/s)');
xlabel('X (m)');
legend('FTCS/FTCS ','Upwind/FTCS','Exact');
title(['Analytical vs Finite Volume Schemes at T=',num2str(T),'s , \Deltax=',num2str(dx),'m , Rec=',num2str(Rec)])
