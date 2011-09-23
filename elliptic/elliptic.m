%% Elliptic PDE

%% 
clc; clear all; close all;

%% Definitions

L = 0.01;			% Length of Plate (x-dir), L (m)
h = 0.02;                       % Height of Plate (y-dir), h (m)
To = 100;                       % Initial Temperature, To (K)
dx = 0.0005;                    % Mesh Size in X-Dir, dx (m)
dy = 0.0005;                    % Mesh Size in Y-Dir, dy (m)
x = [0:dx:L];                   % Create X Vector, x
y = [0:dy:h];                   % Create Y Vector, t
T = zeros(length(x),length(y)); % Create Velocity Matrix

T(:,1) = To;                    % Apply BC: T(x,0)=To
T(1,:) = 0;                     % Apply BC: T(0,y)=0
T(:,end) = 0;                   % Apply BC: T(x,h)=0

beta = dx/dy;                   % Beta, beta
Pomega = 1.191;                 % PSOR Relaxation Term, Pomega
Lomega = 1.29;                  % LSOR Relaxation Term, Lomega
tol = 0.000001;                 % Desired Residual Tolerance, tol


%% Analytical Solution

Ta = zeros(length(x),length(y));
tic;
for i=1:length(x);
    for j=1:length(y);
        SUM=0;
        temp=0;
        for n=1:10000;
            temp=((2*To)/(n*pi*(1-exp(-n*pi*h/L)))) * (1-cos(n*pi)) * (exp(-n*pi*y(j)/(2*L))-exp((-n*pi/L)*(h-((y(j))/2)))) * sin(n*pi*x(i)/(2*L));
            SUM=SUM+temp;
        end
        Ta(i,j)=SUM;
    end
end
Tatime=toc

%% Point Successive Over-Relaxation (PSOR)

PDiff = 1;                      % Initial Value for PSOR Difference, PDiff
PTold=T;                        % Create Numerical Matrix PTold
PTnew=T;                        % Create Numerical Matrix PTnew
PIter=1;                        % Start Iterations at O
Pres=ones(5,1);
tic;
while min(Pres)>tol
   
    for i=2:length(x)-1;
        for j=2:length(y)-1;
          PTnew(i,j)=(1-Pomega)*PTold(i,j)+(Pomega/(2*(1+beta^2)))*(PTold(i+1,j)+PTnew(i-1,j)+ (beta^2)*(PTold(i,j+1)+PTnew(i,j-1)) );
        end
    end
    for j=1:length(y);
          PTnew(end,j)=(4*PTnew(end-1,j)-PTnew(end-2,j))/3;
    end
    Pres(PIter,1)=norm((PTnew(:)-PTold(:)))/norm((PTnew(:)));
    PDiff=max(PTnew(:)-Ta(:));
    
    PTold=PTnew;                % Otherwise continue
    PIter=PIter+1;              % Increment PIter
    
end
PSORtime=toc

%% Line Successive Over-Relaxation (LSOR)

LDiff = 1;                      % Initial Value for LSOR Difference, LDiff
LTold=T;                        % Create Numerical Matrix LTold
LTnew=T;                        % Create Numerical Matrix LTnew
LIter=1;                        % Start Iterations at O
Lres=ones(5,1);
tic;
while min(Lres)>tol

for j=2:length(y)-1;

    for i=1:length(x);
        if i==1
        a(i)=0;
        d(i)=1;
        b(i)=0;
        c(i)=0;
        
        elseif i==length(x);
        a(i)=Lomega;
        d(i)=a(i-1) - 3*b(i-1);
        b(i)=d(i-1) + 4*b(i-1);
        c(i)=c(i-1);
        
        else
        a(i)=Lomega;
        d(i)=-2*(1+beta^2);
        b(i)=Lomega;
        c(i)=-(1-Lomega)*(2*(1+beta^2))*LTold(i,j) - (Lomega*beta^2)*(LTold(i,j+1)+LTnew(i,j-1));
        end
        
    end
    Xp=diag(d)+diag(a(1:end-1),1)+diag(b(2:end),-1);
    yp=c';
    LTnew(:,j)=Xp\yp;
end 
    Lres(LIter,1)=norm((LTnew(:)-LTold(:)))/norm(LTnew(:));
    LDiff=max(LTnew(:)-Ta(:));      % Will break loop if tolerance is met
    LTold=LTnew;                    % Otherwise continue
    LIter=LIter+1;                  % Increment PIter
end

LSORtime=toc

%% Alternating Direction Iteration (ADI)

ADiff = 1;                      % Initial Value for ADI Difference, ADiff
ATold=T;                        % Create Numerical Matrix ATold
ATnew=T;                        % Create Numerical Matrix ATnew
AIter=1;                        % Start Iterations at O
Ares=ones(5,1);

tic;
while min(Ares)>tol

% X-Sweep
for j=2:length(y)-1;

    for i=1:length(x);
        if i==1
        ax(i)=0;
        dx(i)=1;
        bx(i)=0;
        cx(i)=0;
        
        elseif i==length(x);
        ax(i)=1;
        dx(i)=ax(i-1) - 3*bx(i-1);
        bx(i)=dx(i-1) + 4*bx(i-1);
        cx(i)=cx(i-1);
        
        else
        ax(i)=1;
        dx(i)=-2*(1+beta^2);
        bx(i)=1;
        cx(i)=(-beta^2)*(ATold(i,j+1)+ATnew(i,j-1));
        end
        
    end
    Xx=diag(dx)+diag(ax(1:end-1),1)+diag(bx(2:end),-1);
    yx=cx';
    ATnew(:,j)=Xx\yx;
end 

% Y-Sweep
for i=2:length(x)-1;

    for j=1:length(y);
        if j==1
        ay(j)=0;
        dy(j)=1;
        by(j)=0;
        cy(j)=To;
        
        elseif j==length(y);
        ay(j)=0;
        dy(j)=1;
        by(j)=0;
        cy(j)=0;
        
        else
        ay(j)=beta^2;
        dy(j)=-2*(1+beta^2);
        by(j)=beta^2;
        cy(j)=-ATnew(i+1,j)-ATnew(i-1,j);
        end
    end
    
    Xy=diag(dy)+diag(ay(1:end-1),1)+diag(by(2:end),-1);
    yy=cy';
    ATnew(i,:)=Xy\yy;
end 
    Ares(AIter,1)=norm(ATnew(:)-ATold(:))/norm(ATnew(:));
    ADiff=max(ATnew(:)-Ta(:));      % Will break loop if tolerance is met
    ATold=ATnew;                    % Otherwise continue
    AIter=AIter+1;                  % Increment AIter
end

ADItime=toc

%% Plot Results
figure;
hold on
plot(x,PTnew(:,11),'ro--');
plot(x,LTnew(:,11),'m+-');
plot(x,ATnew(:,11),'g*--');
plot(x,Ta(:,11),'-');
ylabel('T (K)');
xlabel('X (m)');
legend('PSOR','LSOR','ADI','Exact');
title(['Analytical vs Finite Difference Schemes at Y=',num2str(y(11)),'m'])

figure;
hold on
plot(x,PTnew(:,21),'ro--');
plot(x,LTnew(:,21),'m+-');
plot(x,ATnew(:,21),'g*--');
plot(x,Ta(:,21),'-');
ylabel('T (K)');
xlabel('X (m)');
legend('PSOR','LSOR','ADI','Exact');
title(['Analytical vs Finite Difference Schemes at Y=',num2str(y(21)),'m'])

figure;
hold on
plot(x,PTnew(:,31),'ro--');
plot(x,LTnew(:,31),'m+-');
plot(x,ATnew(:,31),'g*--');
plot(x,Ta(:,31),'-');
ylabel('T (K)');
xlabel('X (m)');
legend('PSOR','LSOR','ADI','Exact');
title(['Analytical vs Finite Difference Schemes at Y=',num2str(y(31)),'m'])