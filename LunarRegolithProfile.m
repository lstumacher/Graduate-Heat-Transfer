function [wall] = MEGN571_HW3_F2024_starter()
%
%%  MEGN571 Module 3 Assignment Fall 2024:
%        Transient heat conduction: analytical and numerical methods by Prof. Jackson
%        See Module 3 assignment write-up for the problem description
%
close all;
%
%% Set properties of the sintered regolith and constants
wall.k = 0.8;           % constant sintered regolith thermal conductivity [W/m*K]
wall.rho = 2500;       % sintered regolith  density [kg/m^3]
wall.c_p = 650;         % sintered regolith specific heat capacity [J/kg*K]
wall.c_p_var = [325, 1.3];  % sintered regolith specific heat capacity for c_p = c_p(1) + c_p(2)*T with T in K
wall.eps = 0.8;    % wall emissivity
wall.alpha = wall.k/(wall.rho*wall.c_p); % plate thermal diffusivity [m^2/s]
sigma = 5.67e-8;  % Stefan-Boltzman constant [W m^-2 K^-4]
%
%% Set boundary conditions and initial condiions
wall.T_inf  = 90;     %  external temperature to which the wall is radiating [K]
wall.T_HX = 580;    % internal HX temperature which is radiating to the wall [K]
wall.T_init = wall.T_inf;    % initial temperature of the wall [K]
%
%%  Set wall geometry for both solid wall and hollow wall
wall.L_x = 0.3;       % sintered regolith wall thickness [m]
wall.dx = 0.01;  % distance between nodes [m]
wall.x = (0:wall.dx:wall.L_x)';      % sintered regolith x-locations from surface [m] 
wall.x_star = wall.x/wall.L_x;   % plate cell locations dimensionless 
wall.dx_cell = [0.5*wall.dx; wall.dx*ones((length(wall.x)-2),1); 0.5*wall.dx];  % cell lengths [m]
%
wall.x_b = 0.1;        % sintered regolith with hollow core base thickness [m]
wall.x_f  = 0.2;        % sintered regolith x location at end of gap [m]
wall.dy_g = 0.22;        % sintered regolith gap height for hollow core wall [m]
wall.dy_f = 0.08;        % sintered regolith fin height for hollow core wall [m]
wall.x = (0:0.01:wall.L_x)';      % sintered regolith x-locations from surface [m] 
wall.x_star = wall.x/wall.L_x;   % plate cell locations dimensionless
%
%  For solid wall cells, set up conduction areas and cell volumes per unit height
wall.A_cell_solid = (wall.dy_g+wall.dy_f)*ones(length(wall.x)+1,1);  % area per unit height of cell for solid wall [m^2/m]
wall.V_cell_solid = (wall.dy_g+wall.dy_f)*wall.dx_cell; % volume per unit heigh of cell for solid wall [m^3/m]
%
%  For hollow wall cells, set up conduction areas and cell volumes per unit height
i_x = find(wall.x < wall.x_b);
wall.A_cell_hollow(i_x) = wall.dy_g+wall.dy_f;  % area per unit height of cell for solid wall [m^2/m]
wall.V_cell_hollow(i_x) = (wall.dy_g+wall.dy_f)*wall.dx_cell(i_x); % volume per unit height of cell for solid wall [m^3/m]
i_x = find(wall.x == wall.x_b);
wall.A_cell_hollow(i_x) = wall.dy_g+wall.dy_f;  % area per unit height of cell for solid wall [m^2/m]
wall.V_cell_hollow(i_x) = (0.5*wall.dy_g+wall.dy_f)*wall.dx_cell(i_x); % volume per unit height of cell for solid wall [m^3/m]
i_x = find(wall.x > wall.x_b & wall.x < wall.x_f);
wall.A_cell_hollow(i_x) = wall.dy_f;  % area per unit height of cell for solid wall [m^2/m]
wall.V_cell_hollow(i_x) = wall.dy_f*wall.dx_cell(i_x); % volume per unit height of cell for solid wall [m^3/m]
i_x = find(wall.x == wall.x_f);
wall.A_cell_hollow(i_x) = wall.dy_g+wall.dy_f;  % area per unit height of cell for solid wall [m^2/m]
wall.V_cell_hollow(i_x) = (0.5*wall.dy_g+wall.dy_f)*wall.dx_cell(i_x); % volume per unit height of cell for solid wall [m^3/m]
i_x = find(wall.x > wall.x_f);
wall.A_cell_hollow(i_x) = wall.dy_g+wall.dy_f;  % area per unit height of cell for solid wall [m^2/m]
wall.A_cell_hollow(length(wall.x)+1) = wall.dy_g+wall.dy_f;  % area per unit height of cell for solid wall [m^2/m]
wall.V_cell_hollow(i_x) = (wall.dy_g+wall.dy_f)*wall.dx_cell(i_x); % volume per unit height of cell for solid wall [m^3/m]
%
%% Calculate linearized radiative heat transfer coefficient
T_0_guess = wall.T_HX+5; %Guess temperature of left hand side similar to heat exchanger temp
wall.h_rad_HX = sigma*wall.eps*(wall.T_HX^4 - T_0_guess^4)/(wall.T_HX-T_0_guess); %Using initial guess to predict h_rad_HX
T_L_guess = wall.T_inf+5; %Outside edge guess temperature is close to environmental temp
wall.h_rad_inf = sigma*wall.eps*(wall.T_inf^4 - T_L_guess^4)/(wall.T_inf-T_L_guess); %Using rhs temperature guess to predict h_rad_inf
%
%% Calculate the Biot number for the two radiative heat transfers
wall.Biot_HX = wall.h_rad_HX*wall.L_x/wall.k; %Biot Number for HX
wall.Biot_inf = wall.h_rad_inf*wall.L_x/wall.k; %Biot Number for environment
%
%% Set up the times and their dimensionless time (i.e., Fourier number)
dt_step = 60;     %  time interval [s]
t_max = 24*3600;   % maximum time [s]
time = [0:dt_step:t_max];   %  time [s]
wall.time = time; %Establish time as a variable inside the struc
wall.Fo = wall.alpha*wall.time/(wall.L_x^2);  % non-dimensional time or t_star [--]
%
%% Get indices for every 4 hours for plotting
i_t_fig = [0:14400:t_max]/dt_step + 1; %Iterative indices for every 4 hours
%
%% Part a) ---------------------------------------------------------------------
%%  Call the function to find the separation of variable solution with j_max term
j_max = 15; %Number of terms defines calculated in the solutioos
[wall.Temp_exact, wall.qflux_0_exact, wall.qflux_L_exact] = TransientPlate_SepVar(wall, time, j_max); 
%
%% Plot separation of variables temperature profile for 2 hour intervals
figure
plot(wall.x,wall.Temp_exact(:,i_t_fig),'Linewidth',2);
xlabel('Horizontal Position (m)','FontSize',14);
ylabel('Temperature (K)','FontSize',14);
ylim([wall.T_inf wall.T_HX]);
legend(strcat('t = ',num2str(time(i_t_fig(1))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(2))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(3))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(4))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(5))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(6))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(7))/3600),' h')); 
title('Analytical Temperature Profile vs. time','FontSize',14);
%
%-- Part b) -----------------------------------------------------------
%%  Call the function to find the semi-infinite solution
[Temp, wall.qflux_0_seminf] = TransientPlate_Semi_Infinite(wall, time);
%
%%  Find when semi-infinite solution predicts more than 1% difference than separation of variables
i_t_seminf_good = time/dt_step + 1; %Creates indexes in terms of minutes

%
%  Semi-infinite predicts within 1% of heat flux up to time index 652
relative_diff = abs(wall.qflux_0_exact - wall.qflux_0_seminf) ./ abs(wall.qflux_0_exact); %Finds percent difference between solns
diff_indices = find(relative_diff > 0.01); %Isolates indexes where percent diff is greater than 0.01
wall.t_seminf_good = i_t_seminf_good(diff_indices(3)); %Provides first significant index where solutions diverge 
%%  Plot heat flux at x = 0 and x = L_x as a function of time
figure;
% plot(wall.time,wall.qflux_0_seminf(1,:),'g','LineWidth',2)
plot(wall.time,wall.qflux_0_exact(1,:),'r', ...
       wall.time,wall.qflux_0_seminf(1,:),'g',...
       wall.time,wall.qflux_L_exact(1,:),'-r','LineWidth',2)
xlabel('Time [s]','FontSize',14);
ylabel('Wall Heat Flux [W/m^2]','FontSize',14);
ylim([0 5000])
legend('Separation of Variables x = 0', 'Semi-infinite x = 0',...
           'Separation of Variables x = L_x','FontSize',12);
title('Internal and External Wall Heat Flux vs. Time','FontSize',14);
%
%-- Part c -----------------------------------------------------------
%%  Set up the options for the ode solver and the time span fo
options_ode = odeset('RelTol',1e-8,'AbsTol',1e-9,'MaxStep',2400); 
t_span = wall.time;
%
%% Set up the time derivative function handle and the initial temperature
Temp_0 = ones(length(wall.x),1)*wall.T_init; %Guess values for 1D wall
Plate_num = @(t, T)TransientPlate_num(t, wall, T, wall.A_cell_solid, wall.V_cell_solid); %Uses function to set up temperature and time array
%
%%  Run ODE integrator to get wall temperature variation with time
[wall.time_num, wall.Temp_num] = ode45(Plate_num, t_span, Temp_0, options_ode);  %Uses advanced non-linear solver to extract appropriate T vals
wall.Temp_num_solid = wall.Temp_num';
wall.qflux_0_num_solid = wall.h_rad_HX*(wall.T_HX-wall.Temp_num_solid(1,:));
wall.qflux_L_num_solid = wall.h_rad_inf*(-wall.T_inf+wall.Temp_num_solid(end,:)); %Represents heat fluxes across boundaries

%% Plot numerical solution temperature profile for solid wall
figure
plot(wall.x, wall.Temp_num_solid(:,i_t_fig),'Linewidth',2);
xlabel('Horizontal Position (m)','FontSize',14);
ylabel('Temperature (K)','FontSize',14);
ylim([wall.T_inf wall.T_HX]);
legend(strcat('t = ',num2str(time(i_t_fig(1))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(2))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(3))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(4))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(5))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(6))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(7))/3600),' h')); 
title('Numerical Temperature Profile for Solid Wall vs. time','FontSize',14);
%
%%  Plot heat flux at x = 0 and x = L as a function of time for solid wall
figure;
plot(wall.time,wall.qflux_0_num_solid(1,:),'r', ...
       wall.time,wall.qflux_L_num_solid(1,:),'-r','LineWidth',2)
xlabel('Time [s]','FontSize',14);
ylabel('Wall Heat Flux [W/m^2]','FontSize',14);
ylim([0 5000])
legend('Numerical Solid Wall x = 0', 'Numerical Solid Wall x = L_x','FontSize',12);
title('Internal and External Wall Heat Flux vs. Time','FontSize',14);
%
%% Part_d) -----------------------------------------------------------
%% Set up the time derivative function handle and the initial temperature
Temp_0 = ones(length(wall.x),1)*wall.T_init; %Initial Guess
Plate_num = @(t, T)TransientPlate_num(t, wall, T, wall.A_cell_hollow, wall.V_cell_hollow);
%
%%  Run ODE integrator to get wall temperature variation with time
[wall.time_num, wall.Temp_num] = ode45(Plate_num, t_span, Temp_0, options_ode);  
wall.Temp_num_hollow = wall.Temp_num';
wall.qflux_0_num_hollow = wall.h_rad_HX*(wall.T_HX-wall.Temp_num_hollow(1,:));
wall.qflux_L_num_hollow = wall.h_rad_inf*(-wall.T_inf+wall.Temp_num_hollow(end,:));
%
%% Plot numerical solution temperature profile for hollow wall with gap
figure
plot(wall.x, wall.Temp_num_hollow(:,i_t_fig),'Linewidth',2);
xlabel('Horizontal Position (m)','FontSize',14);
ylabel('Temperature (K)','FontSize',14);
ylim([wall.T_inf wall.T_HX]);
legend(strcat('t = ',num2str(time(i_t_fig(1))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(2))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(3))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(4))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(5))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(6))/3600),' h'), ...
           strcat('t = ',num2str(time(i_t_fig(7))/3600),' h')); 
title('Numerical Temperature Profile for Wall with Gap vs. time','FontSize',14);
%
%%  Plot heat flux at x = 0 and x = L as a function of time for hollow wall with gap
figure;
plot(wall.time,wall.qflux_0_num_hollow(1,:),'r', ...
       wall.time,wall.qflux_L_num_hollow(1,:),'-r','LineWidth',2)
xlabel('Time [s]','FontSize',14);
ylabel('Wall Heat Flux [W/m^2]','FontSize',14);
ylim([0 5000])
legend('Numerical Wall with Gap x = 0', 'Numerical wall with Gap x = L_x','FontSize',12);
title('Internal and External Wall Heat Flux vs. Time','FontSize',14);
end
%
%--------------------------------------------------------
%%  Function to calculate separation of variable solution
function [Temp, qflux_0, qflux_L] = TransientPlate_SepVar(wall, time, j_max)
%
%  Allocate memory for dimensionless temeprature and its derivative at x = 0
theta = zeros(length(wall.x),length(time));
theta_h = zeros(length(wall.x),length(time));
dtheta_dx_star_0 = zeros(1,length(time));
dtheta_dx_star_L = zeros(1,length(time));
%
% Get dimensionless time
Fo = wall.Fo;
%
%  Determine the particular solution which is independent of t
%     because the boundary conditions are constant.

A = [wall.Biot_HX, -1; wall.Biot_inf, 1 + wall.Biot_inf]; % Uses boundary conditions to create a matrix for LHS of cp coefficients
b = [wall.Biot_HX; 0]; % Uses boundary conditions to create RHS vector that equals A*cp 
cp = A\b; % Solves for cp
theta_p = cp(1) + cp(2)*wall.x_star; % Solution for guess value of theta particular

% Solve for the eigenvalues lambda_n of the 1-D transient conduction homogeneous solution  (Part_a) + 0.4 pts.)
%    with mixed boundary conditions at x_star = 0 and x_star = 1.   
lambda_n = zeros(j_max,1);
for j = 1:j_max %Iterates through number of terms
    lambda_0 = pi*(j-1/2); % Defines guess values for lambda given where tan(lx) tends to zero
    func = @(x) (-x+wall.Biot_HX*wall.Biot_inf/x)*sin(x)+(wall.Biot_HX+wall.Biot_inf)*cos(x);
    lambda_n(j) = fzero(func,lambda_0); %Sets up and solves for lambda in the equation derived from homogenous BC
end
%
%  Find initial conditions for homogeneous solution to get constant c_j 
for j = 1:j_max %Iterates through number of terms
    funcnum = @(x,j) (-cp(1) -cp(2)*x) .* (cos(lambda_n(j).*x)+wall.Biot_HX*sin(lambda_n(j).*x)/lambda_n(j)); %Num of IC/orthogonality
    funcden = @(x,j) (cos(lambda_n(j).*x)+wall.Biot_HX.*sin(lambda_n(j).*x)/lambda_n(j)).^2; %Denominator of IC/orthogonality eq
    num = integral(@(x) funcnum(x,j), 0,1); 
    bum = integral(@(x) funcden(x,j), 0, 1);
    C_j(j) = num/bum; % Solves for C_1j constant
end
%
% Calculate the sum for the  transient 1-D solution from separation of variables 

theta_h = zeros(length(wall.x),length(time)); %Predefines theta_h 

for x_star = 1:length(wall.x_star) %Solves for theta h at every time and positon using constants derived from BCs and IC.
     for t = 1:length(Fo)
         for j = 1:j_max
            theta_h(x_star,t) =theta_h(x_star,t)+C_j(j)*exp(-lambda_n(j).^2.*Fo(t)).*(cos(lambda_n(j).*wall.x_star(x_star))+(wall.Biot_HX/lambda_n(j)).*sin(lambda_n(j).*wall.x_star(x_star))); 
         end
     end
end

 theta = theta_h + theta_p; %Adding solutions together
%
%  Translate theta to temperature and dtheta_stardx_star to heat flux
Temp =  theta * (wall.T_HX - wall.T_inf) + wall.T_inf; % Finds T from selected theta equation
qflux_0 =  wall.h_rad_HX*(wall.T_HX-Temp(1,:)); % Models the radiative heat coming in through the LHS
qflux_L = wall.h_rad_inf*(-wall.T_HX+Temp(end,:)); % Models radiative heat leaving to environment 
end
%
%--------------------------------------------------------
%%  Function to calculate semi-infinite plate solution
function [Temp, qflux_0_seminf] = TransientPlate_Semi_Infinite(wall, time)
%
% Solve the semi-infinite plate solution
%  and calculate the surface and midplane values of temperature.
%  Here reverse x to make it consistent with x for 1-D transient conduction.
%
% Allocate memory for dimensionless temeprature and its derivative at x = 0
theta_star = zeros(length(wall.x),length(time));
dtheta_stardx_0 = zeros(1,length(time));
%
% 
Fo = wall.Fo; % Non dimensionalized time
 % 

  for x_star = 1:length(wall.x) %Uses provided semi-infite equation to find temperature at all points up to the midplane in time
    for t = 1:length(time)
            theta_star(x_star, t) =  theta_star(x_star, t) + erfc(wall.x(x_star)./(2*sqrt(wall.alpha*time(t)))) -...
                exp(wall.h_rad_HX*wall.x(x_star)/wall.k + (wall.h_rad_HX^2)*(wall.alpha*time(t))/(wall.k^2)).*erfc(wall.x(x_star)./(2*sqrt(wall.alpha*time(t))) + (wall.h_rad_HX/wall.k)*sqrt(wall.alpha*time(t)));
     end
 end

%
%  Translate theta_star to temperature and dtheta_stardx_star to heat flux (Part_b) + 0.1 pts.)
Temp = theta_star*(wall.T_HX-wall.T_inf)+wall.T_inf;
qflux_0_seminf = wall.h_rad_HX*(wall.T_HX-Temp(1,:)); %Head entering through x=0 via radiative flux
end
%
%--------------------------------------------------------
%% Function for calculating numerical time derivatives for wall temperature profile
function [dTdt] = TransientPlate_num(t, wall, T,A,V)
%
% Allocate meory for temperature derivative
dTdt = zeros(length(wall.x),1);
%
%  Get variable c_p value
c_p = wall.c_p_var(1) + wall.c_p_var(2)*T; %Data shows c_p of regolith varies linearly with temp
%
%
for i = 2:length(wall.x)-1 %Iterating through central points in 1D
    dTdt(i) = (wall.k*A(i+1)*(T(i+1)-T(i))/wall.dx + wall.k*A(i)*(T(i-1)-T(i))/wall.dx) / (V(i)*wall.rho * c_p(i)); %Finite volume approximation
end
dTdt(1) = (wall.k*A(2)*(T(2)-T(1))/wall.dx + wall.h_rad_HX * A(1)*(wall.T_HX - T(1)))/(V(1)*wall.rho*c_p(1)); %BC for heat entering LHS
dTdt(end) = (wall.k*A(end-1)*(T(end-1)-T(end))/wall.dx + wall.h_rad_inf * A(end)*(-wall.T_inf + T(end)))/(V(end)*wall.rho*c_p(end)); %BC for heating leaving RHS
%
end
