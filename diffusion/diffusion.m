%% Diffusion PDE

%% 
clc; clear all; close all;

%% Definitions

h = 0.04;	% Distance, h (m)
T = 1.08;	% Time, T (s)
nu = 0.000217;	% Kinematic Viscosity, nu (m^2/s)
rho = 800;	% Density, rho (kg/m^3)
Uo = 40;	% Speed of Lower Plate, Uo (m/s)

%% Computational Method

dy = 0.001;	% Step Size in Y, dy

% Step Size in T (Calculated), dt
dt = (0.5*dy^2)/nu;	

% Constant of Propogation, d
d = nu*dt/dy^2;		

% Create Y Column Vector, y
y = [0:dy:h]';

% Create Time Step Row Vector, t
t = [0:dt:T];

% Number of Time Steps, numT
numT = length(t);

% Number of Y-Direction Mesh Points, numY
numY = length(y);

% Create Velocity Matrix
U = zeros(numY,numT);

% Apply Initial Condition
U(1,:) = Uo;

% Solve
for j=2:numT
    for i=2:numY-1
        U(i,j) = d*U(i-1,j-1) + (1-2*d)*U(i,j-1) + d*U(i+1,j-1);
    end
end

%% Analytical Solution
ya = 0:0.001:h;

eta = ya/(2*sqrt(nu*T));
eta1 = h/(2*sqrt(nu*T));
 
SUM1 = 0;
SUM2 = 0;

for n=0:10000
    X1=erfc(2*n*eta1+eta);
    SUM1=SUM1+X1;
end

for n=1:9999
    X2=erfc(2*n*eta1-eta);
    SUM2=SUM2+X2;
end

Ua = Uo*(SUM1-SUM2)';

%% Plot Results
figure;
plot(U(:,numT),y,'ro--');
hold on
plot(Ua,ya,'-');
xlabel('U (m/s)');
ylabel('Y (m)');
legend('FTCS ','Exact');
title(['Analytical vs FTCS Scheme at T=',num2str(T),'s , \Deltay=',num2str(dy),'m, \Deltat=',num2str(dt),'s' ])

