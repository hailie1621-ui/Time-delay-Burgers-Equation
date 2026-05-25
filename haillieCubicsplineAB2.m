%========================================================================================================
% project name- Numerical Solution for Time-Delay Burgers' Equation Using Combined Cubic spline AB(4nd) Method
%========================================================================================================
% By: HAILEYESUS TIGABIE ASMARE,  ID Number DTU14Ms023
%========================================================================================================
% DEBRETABOR UNIVERSITY COLLEGE OF NATURAL AND COMPUTATIONAL SCIENCE DEPARTMENT OF MATHEMATICS
%========================================================================================================
% NUMERICAL ANALYSIS Stream For Masters Summer Program 2018 EC
%========================================================================================================
tic
clc; clear all; close all;
L = input('Enter the interval of space  L = ');   % let x in [0, L]
T = input('Enter the interval of time  T = ');    % let t in [0, T]
h = input('Enter the step size in space  h = ');  % 
h1 = h/2;
k = input('Enter the step size in time  k = ');   % apply CLF stability condition 
z = input('Enter the time delay (tau) z = '); 
v = input('Enter the viscosity coefficient  v = ');
Nx = round(L/h);            % number of grids in space
Nx1 = round(L/h1);          % second season number of grids in space
Nt = round(T/k);            % number of grids in time for both season
r = z/k;                    % terms of delay for both season
b=Nt;                       % the last grid for calculate maximum error
ds = round(r);              % delay index 
% GRID
x = linspace(0,L,Nx+1);      % 
t = linspace(0,T,Nt+1);
x1 = linspace(0,L,Nx1+1);    
% ------------------------------
% INITIAL CONDITION
% ------------------------------
u = zeros(Nx+1, Nt+1);   
for i = 1:Nx+1
    u(i,1) = sin(pi*x(i)); %x(i); %(2*v*pi*sin(pi*x(i)))/(2+cos(pi*x(i))); history function 
end
u1 = zeros(Nx1+1, Nt+1);   
for i = 1:Nx1+1
    u1(i,1) = sin(pi*x1(i)); %x1(i); %(2*v*pi*sin(pi*x1(i)))/(2+cos(pi*x1(i)));% history function 
end
% ------------------------------
% BOUNDARY CONDITION (Dirichlet Boundary)
% ------------------------------
for j = 1:Nt+1
    u(1,j) = 0;        % u(0,t)
    u(Nx+1,j) = 0; % 1/(1+t(j));         
end
for j = 1:Nt+1
    u1(1,j) = 0;        % u(0,t)
    u1(Nx1+1,j) = 0;  % 1/(1+t(j));  
end
% ------------------------------
% CUBIC SPLINE MATRIX
% ------------------------------
B = zeros(Nx+1,Nx+1);      % coefficient matrices B for first season
rhs = zeros(Nx+1,1);       % rhs is right hand side for first season
A = zeros(Nx1+1,Nx1+1);     % coefficient matrices A for second season 
rhs1 = zeros(Nx1+1,1);     % rhs is right hand side for second season
B(1,1) = 1;               
B(Nx+1,Nx+1) = 1;          
for i = 2:Nx
    B(i,i-1) = h;          
    B(i,i)   = 4*h;        
    B(i,i+1) = h;           
end
A(1,1) = 1;               
A(Nx1+1,Nx1+1) = 1;          
for i = 2:Nx1
    A(i,i-1) = h1;          
    A(i,i)   = 4*h1;        
    A(i,i+1) = h1;           
end
% ------------------------------
% TIME STEPPING
% ------------------------------
f = zeros(Nx+1, Nt+1);      
f1 = zeros(Nx1+1, Nt+1); 
 for j = 1:Nt
    % --- RHS for spline second derivative ---
    rhs = zeros(Nx+1,1);
    for i = 2:Nx
        rhs(i) = 6*(u(i+1,j) - 2*u(i,j) + u(i-1,j))/h;
    end
    M = B\rhs;             %(M = S'')
    % --- Time delay handling ---
    if abs(r-ds) <= 1.0e-10    
      m = ds;
      if j <= m
        ud = u(:,1);     % History function (t <= tau  u(x,0)
      else
        ud = u(:,j-m);
      end
    else           
       m = floor(r);
       theta = r - m;
      if j <= m+1
        ud = u(:,1);
      else
         ud = (1-theta) * u(:,j-m) + theta * u(:,j-m-1);
      end
    end
    % --- First step (RK-2) ---
    if j == 1
        k1 = zeros(Nx+1,1);
        for i = 2:Nx
            ux = (u(i+1,j) - u(i-1,j))/(2*h) - (h/6)*(M(i+1) - M(i-1));
            k1(i) = -ud(i)*ux + v*M(i);
        end
        ut = u(:,j); 
        for i = 2:Nx
            ut(i) = u(i,j) + k*k1(i);
        end
        % ---- Spline for temporary solution ----
        rhs_t = zeros(Nx+1,1);
        for i = 2:Nx
            rhs_t(i) = 6*(ut(i+1) - 2*ut(i) + ut(i-1))/h;
        end
        Mt = B\rhs_t;     
        % ---- k2 ----
        k2 = zeros(Nx+1,1);
        for i = 2:Nx
            utx = (ut(i+1) - ut(i-1))/(2*h) - (h/6)*(Mt(i+1) - Mt(i-1));
            k2(i) = -ud(i)*utx + v*Mt(i);
        end
        % Final update for j=1 ----
        for i = 2:Nx
            u(i,j+1) = u(i,j) + (k/2)*(k1(i) + k2(i));
        end
    % --- Adams-Bashforth (2-step) ---
    else
        for i = 2:Nx
            ux = (u(i+1,j) - u(i-1,j))/(2*h) - (h/6)*(M(i+1) - M(i-1));
            f(i,j) = -ud(i)*ux + v*M(i);
            
            % 2nd order Adams-Bashforth ????
            u(i,j+1) = u(i,j) + (k/2)*(3*f(i,j) - f(i,j-1));
        end
    end
    
 end
% ====================== %
% SECOND GRID (h1 = h/2) %
% ====================== %
 for j=1:Nt
    % --- RHS for spline second derivative ---
    rhs1 = zeros(Nx1+1,1);
    for i = 2:Nx1
        rhs1(i) = 6*(u1(i+1,j) - 2*u1(i,j) + u1(i-1,j))/(h1);
    end
    M1 = A\rhs1;             %  (M = S'')
    % --- Time delay handling ---
    if abs(r-ds) <= 1.0e-10    
      m = ds;
      if j <= m
        ud1 = u1(:,1);     % History function (t <= tau  u(x,0)  )
      else
        ud1 = u1(:,j-m);
      end
    else           
       m = floor(r);
       theta = r - m;
      if j <= m+1
        ud1 = u1(:,1);
      else
         ud1 = (1-theta) * u1(:,j-m) + theta * u1(:,j-m-1);
      end
    end
    % --- First step (RK-2) ---
    if j == 1
        k1 = zeros(Nx1+1,1);
        for i = 2:Nx1
            u1x = (u1(i+1,j) - u1(i-1,j))/(2*h1) - (h1/6)*(M1(i+1) - M1(i-1));
            k1(i) = -ud1(i)*u1x + v*M1(i);
        end
        u1t = u1(:,j); 
        for i = 2:Nx1
            u1t(i) = u1(i,j) + k*k1(i);
        end
        % ---- Spline for temporary solution ----
        rhs_t1 = zeros(Nx1+1,1); % temporary right hand side
        for i = 2:Nx1
            rhs_t1(i) = 6*(u1t(i+1) - 2*u1t(i) + u1t(i-1))/h1;
        end
          M1t = A\rhs_t1;     % temporary splin martices 
          % ---- k2 ----
          k2 = zeros(Nx1+1,1);
        for i = 2:Nx1
            u1tx = (u1t(i+1) - u1t(i-1))/(2*h1) - (h1/6)*(M1t(i+1) - M1t(i-1));
            k2(i) = -ud1(i)*u1tx + v*M1t(i);
        end
            % Final update for j=1 ----
        for i = 2:Nx1
            u1(i,j+1) = u1(i,j) + (k/2)*(k1(i) + k2(i));
        end
    % --- Adams-Bashforth (2-step) ---
    else
        for i = 2:Nx1
            u1x = (u1(i+1,j) - u1(i-1,j))/(2*h1) - (h1/6)*(M1(i+1) - M1(i-1));
            f1(i,j) = -ud1(i)*u1x + v*M1(i);
            % 2nd order Adams-Bashforth ????
            u1(i,j+1) = u1(i,j) + (k/2)*(3*f1(i,j) - f1(i,j-1));
        end
    end
 end
  
% maximum error between two mesh size h and h/2
 Error = max(abs(u1(1:2:end, b) - u(:, b)));
fprintf('Absolute maximum Error is\n  %1.6e\n', Error);

% ------------------------------
% GRAPHICAL REPRESENTATION
% ------------------------------
figure; hold on;
plot(x1, u1(:,1), 'k-', 'LineWidth', 2.5);                     % t = 0
plot(x1, u1(:, round(Nt/8)), 'y-', 'LineWidth', 2.5);
plot(x1, u1(:, round(Nt/4)), 'r-', 'LineWidth', 2.5);         
plot(x1, u1(:, round(Nt/2)), 'b-', 'LineWidth', 2.5);         
plot(x1, u1(:, end), 'g-', 'LineWidth', 2);                  % Final Time

xlabel('x1');
ylabel('u1(x,t)');
title('Solution at Different Time Levels');
legend('t = 0','t=T/8','t = T/4','t = T/2','Final Level (t = T)');
grid on;

figure;
[X,Tm] = meshgrid(x1,t);
surf(X,Tm,u1')
xlabel('x'); ylabel('t'); zlabel('u1(x,t)')
title('Cubic Spline-Adams Bashforth Solution 3D')
shading interp
colorbar;

toc
