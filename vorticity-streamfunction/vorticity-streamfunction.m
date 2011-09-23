%% Vorticity-Stream Function Method

%% 
clc; clear all; close all;

%% Definitions

L = 5.0;			% Length of Channel (x-dir), L
h = 1.0;                        % Height of Channel (y-dir), h
Nx = 50;                        % Number of Grid Points (x-dir), Nx
Ny = 20;                        % Number of Grid Points (y-dir), Ny
dx = L/(Nx-1);                  % Mesh Size in X-Dir, dx
dy = h/(Ny-1);                  % Mesh Size in Y-Dir, dy
x = [0:dx:L];                   % Create X Vector, x
y = [0:dy:h];                   % Create Y Vector, t

Re = 10;                        % Reynolds Number, Re
T = 0;                          % Initial Time, T
dt = 0.001;                     % Time Step, dt
tol = 0.0000001;                % Desired Tolerance, tol
omega = 1.91;                   % Relaxation Term, omega
beta = dx/dy;                   % Mesh Size Relation, beta

U = ones(Nx,Ny);                % Create U - Velocity Matrix
V = zeros(Nx,Ny);               % Create V - Velocity Matrix
Xi = zeros(Nx,Ny);              % Create Vorticity Matrix
Psi = zeros(Nx,Ny);             % Create Stream-Function Matrix

% Inlet
U(1,:) = 1;                     % Apply BC: U(0,j)=1
V(1,:) = 0;                     % Apply BC: V(0,j)=0
Xi(1,:) = 0;                    % Apply BC: Xi(0,j)=0
Psi(1,:) = y;                   % Apply BC: Psi(0,j)=y(j)

% Outlet
U(Nx,:) = U(Nx-1,:);            % Apply BC: U(Nx,j)=U(Nx-1,j)
V(Nx,:) = V(Nx-1,:);            % Apply BC: V(Nx,j)=V(Nx-1,j)
Xi(Nx,:) = Xi(Nx-1,:);          % Apply BC: Xi(Nx,j)=Xi(Nx-1,j)
Psi(Nx,:) = Psi(Nx-1,:);        % Apply BC: Psi(Nx,j)=Psi(Nx-1,j)

% No-Slip Wall (y=0)
Psi(:,1) = 0;                   % Apply BC: Psi(i,0)=0
Xi(:,1) = (7*Psi(:,1)-8*Psi(:,2)+Psi(:,3))/(2*dy^2);

% No-Slip Wall (y=h)
Psi(:,Ny) = 1;                  % Apply BC: Psi(i,Ny)=1
Xi(:,Ny) = (7*Psi(:,Ny)-8*Psi(:,Ny-1)+Psi(:,Ny-2))/(2*dy^2);


%% Analytical Solution

Ua = zeros(Nx,Ny);

for i=1:Nx;
    for j=1:Ny;
        Ua(i,j)=6*(y(j)-y(j)^2);
    end    
end

%% Numerical Approximation

del = 1;
PSORdel = 1;
Xinew=Xi;
Xiold=Xi;
Psinew=Psi;
Psiold=Psi;
Uold=U;
Unew=U;
Vold=V;
Vnew=V;

while del>tol
    
    % FTCS
    for i=2:Nx;
        for j=1:Ny;
            if i==Nx
                Xinew(i,j) = Xinew(i-1,j);
            elseif j==1
                Xinew(i,j) = (7*Psinew(i,j)-8*Psinew(i,j+1)+Psinew(i,j+2))/(2*dy^2);
            elseif j==Ny
                Xinew(i,j) = (7*Psinew(i,j)-8*Psinew(i,j-1)+Psinew(i,j-2))/(2*dy^2);
            else
                Xinew(i,j)=Xiold(i,j)...
                + dt*( (1/Re)*( ((Xiold(i+1,j) - 2*Xiold(i,j) + Xiold(i-1,j))/(dx^2)) ...
                + ((Xiold(i,j+1) - 2*Xiold(i,j) + Xiold(i,j-1))/(dy^2)) )...
                - ((Uold(i+1,j)*Xiold(i+1,j) - (Uold(i-1,j)*Xiold(i+1,j) ))/(2*dx))...
                - ((Vold(i,j+1)*Xiold(i,j+1) - (Vold(i,j-1)*Xiold(i,j-1) ))/(2*dy)) );
            end
        end
    end
    
    % PSOR
    while PSORdel>tol
        for i=2:Nx;
            for j=1:Ny;
                if i==Nx
                    Psinew(i,j)=Psinew(i-1,j);
                elseif j==1
                    Psinew(i,j)=0;
                elseif j==Ny
                    Psinew(i,j)=1;
                else
                    Psinew(i,j)=(1-omega)*Psiold(i,j)+(omega/(2*(1+beta^2)))*...
                    ((dx^2)*Xinew(i,j)+Psiold(i+1,j)+Psinew(i-1,j)+(beta^2)*(Psiold(i,j+1)+Psinew(i,j-1)) );
                end
            end
        end
        PSORdel=norm((Psinew(:)-Psiold(:)))/norm(Psiold(:));
        Psiold=Psinew;
    end
    
    % U Matrix
    for i=2:Nx;
        for j=1:Ny;
            if i==Nx
                Unew(i,j)=Unew(i-1,j);
            elseif j==1
                Unew(i,j)=0;
            elseif j==Ny
                Unew(i,j)=0;
            else
                Unew(i,j)=(Psinew(i,j+1)-Psinew(i,j-1))/(2*dy);
            end
        end
    end
    
    % V Matrix
    for i=2:Nx;
        for j=1:Ny;
            if i==Nx
                Vnew(i,j)=Vnew(i-1,j);
            elseif j==1
                Vnew(i,j)=0;
            elseif j==Ny
                Vnew(i,j)=0;
            else
                Vnew(i,j)=-(Psinew(i+1,j)-Psinew(i-1,j))/(2*dx);
            end
        end
    end

    del=norm((Xinew(:)-Xiold(:)))/norm(Xiold(:));
    Xiold=Xinew;
    Psiold=Psinew;
    Uold=Unew;
    Vold=Vnew;
    T=T+dt;
    PSORdel=1;
end

%% Determine Entrance Length

el=1;
eldel=1;

while eldel>0.001
    eldel=norm((Unew(el+1,:)-Unew(el,:)))/norm(Unew(el,:));
    el=el+1;
end

%% Plot Results
figure;
hold on
plot(Unew(el,:),y,'ro-');
plot(Ua(el,:),y,'-');
ylabel('Y');
xlabel('u-velocity');
legend('Approximation','Exact');
title(['Analytical vs FD Approximation: Fully Developed Flow at x=',num2str(x(el)),' and t=',num2str(T)])

figure;
hold on
plot(Unew(2:el,:),y);
ylabel('Y');
xlabel('u-velocity');
legend(['x=',num2str(x(el-6))],['x=',num2str(x(el-5))],['x=',num2str(x(el-4))],['x=',num2str(x(el-3))],['x=',num2str(x(el-2))],['x=',num2str(x(el-1))],['x=',num2str(x(el))]);
title(['Steady-State u-velocity Profile'])
