function [hay, plastic, air] = MEGN571_HW1_F2024_starter()
%
%% MEGN571 HW #1 Fall 2024 starter code: Hay bale oxidation and fire prevection
% *Problem statement*  for hay bale heating due to decomposition is in MEGN571_HW1_S2024.pdf
close all;   % This line closes figures that are open.  It is not necessary but it can be helpful for reporting.
%
%%  Set the geometric, themodynamic, and transport properties
%     for the hay bale, polyethylene wrapper, and ambient
%
%  Set geometric z and r meshes
hay.d = 1.8;     % bale diameter [m]
hay.z_max =  2.4;     % bale length [m]
plastic.t = 0.002;  % plastic wrap thickness [m]
%
hay.dz = 0.01;   % axial discretization for plotting and for numerical solution [m]
hay.z = (0:hay.dz:0.5*hay.z_max); % hay axial locations for plotting and for numerical solution [m]
hay.n_z = length(hay.z);   % number of nodes in z mesh
%
hay.r_max = 0.5*hay.d; % bale radius [m]
hay.dr = 0.01;   % radial discretization for plotting and for numerical solution [m]
hay.r = (0:hay.dr:hay.r_max);  % hay radial locations for plotting and for numerical solution [m]
hay.n_r = length(hay.r);   % number of nodes in r mesh
%
%  Define the bale thermal properties for the bale  (loosely and densely packed) and the plastic
hay.k_a = [0.02 0.06];   % both hay bale thermal conductivities [W/m*K]
hay.rho_a = [100 180]; % both hay bale densities [kg/m^3]
hay.c_v = 4000;  % hay specific heat capacity [J/kg*K]
plastic.k = 0.50; % plastic wrap thermal conductivity [W/m*K]
plastic.eps = 0.84;  % plastic emissivity [--]
plasticext.T = 303.15; %Temperature of the external surface (K)
%
% Define constants and surrounding temperatures
sigma = 5.67e-8;  % Stefan-Boltzmann constant [W m^-2 K^-4]
air.T = 25 + 273.15;        % ambient air temperature [K]
air.h = 12;     % ambient heat transfer coefficient [W m^-2 K^-1]
%
%  Set the heat generation parameters
hay.q_reac_a = [5.0];   % hay decomposition heat generation rates for part 1  [W m^-2]
hay.T_ign = 350;  % hay ignition temperature [K]
hay.T_a_reac = 4100;  % hay decomposition activation temperature [K]
hay.A_reac = 3.8e6;   %  hay decomposition pre-exponential [W m^-3]
%
%%  *Part a) Solution:*
%  Calculate total thermal resistance for the external boundary condition 

h_rad = (sigma*plastic.eps)*(plasticext.T^2+air.T^2)*(plasticext.T+air.T);  % effective radiation heat transfer coefficient [W m^-2 K^-2]
R_cond_plastic = plastic.t/plastic.k; % plastic conduction area-specificresistance [m^2 K W^-1]
R_conv = 1/air.h;
R_rad = 1/h_rad;
hay.R_t_ext = R_cond_plastic+(h_rad+air.h)^-1;  % external resistance [m^2 K W^-1]
%
%  Allocate memory for two analytical solutions -- both temperature
%  profiles and max q_reac
hay.T_z_a = zeros(length(hay.z),length(hay.k_a));  
hay.q_reac_z_max_a = zeros(1,length(hay.k_a));  
hay.T_r_a = zeros(length(hay.r),length(hay.k_a));  
hay.q_reac_r_max_a = zeros(1,length(hay.k_a));  
%
%  Write a comment showing the analytical T(z) with constant q_reac
%  T(z) = -qdot*z^2/(4k) + C1*z + C2
%   \\T(z) = (hay.q_reac_a*0.5*hay.z_max/(2*hay.k_a))*(1-(z/(0.5*hay.z_max))^2) + (hay.q_reac_a*0.5*hay.z_max/)*hay.R_t_ext + air.T
%   B.C. @ z = 0   dT/dz = 0
%   B.C.  @ z = z_max/2  setting q out equal to qcircuit -kdT/dz = T-Tinf/Rth 
hay.q_reac_a = 5.0;  %  thermal energy generation for part a [W m^-3]
for i_k = 1:length(hay.k_a)
    % ???

    hay.T_z_a(:,i_k) = (hay.q_reac_a*0.5*hay.z_max/(2*hay.k_a(i_k)))*(1-(hay.z/(0.5*hay.z_max)).^2) + (hay.q_reac_a*0.5*hay.z_max)*hay.R_t_ext + air.T;   % steady-state axial temperature profile as a function of z [K]
%
%  Find q_reac to give T_ign at z = 0 for analytical axial solution
    hay.q_reac_z_max_a(:,i_k) = (hay.T_ign-air.T) / (((0.5*hay.z_max)^2 / (2 * hay.k_a(i_k))) + (0.5*hay.z_max / 2) * hay.R_t_ext);  % 
%   
end 

%
%  Write a comment showing the analytical T(r) with constant q_reac
%  T(r) = -qdot*r^2/(4k) + C1*ln(r) + C2
%  T(r) = (hay.q_reac_a*hay.r_max/(4*hay.k_a))*(1-(r/(hay.r_max))^2) + (hay.q_reac_a*hay.r_max/(2))*hay.R_t_ext + air.h;
%   B.C. @ r = 0  Adiabatic Centerline: q -> dT/dr = 0 -> C1 = 0
%   dT/dr at r= 0 is zero
%   B.C. @ r = r_max/2 qout = qcircuit allowing you to substitute in Tzbale
%   = T1 = qdot*d/4*Rth + Tinf
%   
for i_k = 1:length(hay.k_a)
    hay.T_r_a(:,i_k) = (hay.q_reac_a*(hay.r_max^2)/(4*hay.k_a(i_k)))*(1-(hay.r/(hay.r_max)).^2) + (hay.q_reac_a*hay.r_max/(2))*hay.R_t_ext + air.T; % steady-state radial temperature profile as a function of r [K]
%   Find q_reac to give T_ign at r = 0 for analytical radial solution
    hay.q_reac_r_max_a(1,i_k) = (hay.T_ign-air.T) / ((hay.r_max^2 / (4 * hay.k_a(i_k))) + (hay.r_max / 2) * hay.R_t_ext); 

end

%  Plot temperature vs. z for analytical axial solution
figure;
plot(hay.z,hay.T_z_a,'LineWidth',2);  
xlabel('width from center plance [m]','Fontsize',14);
ylabel('Temperature, [K]','Fontsize',14);
legend(strcat('loose k =',num2str(hay.k_a(1))), ...
           strcat('dense k =',num2str(hay.k_a(2))), 'Fontsize',12);  
title('Axial Temperature Profiles with Constant Heat Generation')
%
%  Plot temperature vs. r  for analytical radial solution
figure;
plot(hay.r,hay.T_r_a,'LineWidth',2);     
xlabel('radius [m]','Fontsize',14);
ylabel('Temperature, [K]','Fontsize',14);
legend(strcat('loose k =',num2str(hay.k_a(1))), ...
           strcat('dense k =',num2str(hay.k_a(2))), 'Fontsize',12);  
title('Radial Temperature Profile with Constant Heat Generation')
%
%%  *Part b) Solution:*
%  Set initial guess for two numerical solutions
T_z_0 = hay.T_z_a(:,2); %Setting initial guess equal to steady state conditions
T_r_0 = hay.T_r_a(:,2);
%
% Set up appropriate k and rho for dense hay bale
hay.k_b = hay.k_a(end);
hay.rho_b = hay.rho_a(end);
%
% Set up options for fsolve to find the steady temperature profile
%  Choices for fsolve algorithm -- 'trust-region-dogleg' (default), 'trust-region', and 'levenberg-marquardt'.
% options_solve = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', 'display','iter',...
%    'SpecifyObjectiveGradient' ,false, 'MaxFunEvals', 20000, 'MaxIter', 200, 'TolFun', 1e-8, 'StepTolerance', 1e-12);
%
% Set up options for fsolve to find the steady temperature profile
%  Choices for fsolve algorithm -- 'trust-region-reflective' (default) and 'levenberg-marquardt'.
options_solve = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'display','iter',...
    'SpecifyObjectiveGradient', false, 'MaxFunEvals', 20000,'MaxIter', 200, 'TolFun', 1e-8, 'StepTolerance', 1e-12);
%
% Set up function handles for the axial and radial numerical 
Temp_axial_ss = @(T)Temp_axial(T, hay, air);          % (Part b) +0.2 pts.)
Temp_radial_ss = @(T)Temp_radial(T, hay, air);        % (Part b) +0.2 pts.)
%
% Call fsolve for numerical solution of axial and radial profiles
%   with temperature-depedent heat generation
% hay.T_z_b = fsolve(Temp_axial_ss, T_z_0, options_solve);   % (Part b) +0.2 pts.)
% hay.T_r_b = fsolve(Temp_radial_ss, T_r_0, options_solve);   % (Part b) +0.2 pts.)
%
% Call lqsqnonlin for numerical solution of axial and radial profiles
%   with temperature-depedent heat generation after setting up lower and
%   upper bounds for solutions.
lb_z = air.T*ones(hay.n_z,1); 
ub_z = 2000*ones(hay.n_z,1); 
lb_r = air.T*ones(hay.n_r,1); 
ub_r = 2000*ones(hay.n_r,1); 
%
hay.T_z_b = lsqnonlin(Temp_axial_ss, T_z_0, lb_z, ub_z, options_solve); 
hay.T_r_b = lsqnonlin(Temp_radial_ss, T_r_0, lb_r, ub_r, options_solve);  
%
%  Plot temperature vs. z for numerial axial solution with temperature-depedent heat generation
figure;
plot(hay.z,hay.T_z_b,'LineWidth',2);   
xlabel('width from center plance [m]','Fontsize',14);
ylabel('Temperature, [K]','Fontsize',14);
title('Axial Temperature Profile with T-dependent Heat Generation')
%
%  Plot temperature vs. z for numerial radial solution with temperature-depedent heat generation
figure;
plot(hay.r,hay.T_r_b,'LineWidth',2);   %
xlabel('radius [m]','Fontsize',14);
ylabel('Temperature, [K]','Fontsize',14);
title('Radial Temperature Profile with T-dependent Heat Generation')
%
%%  *Part c) Solution:*
%  Set up initial condition for transient numerical solutions
T_z_t_0 = air.T*ones(length(hay.z),1);    % Initialize with ones
T_r_t_0 = air.T*ones(length(hay.r),1);     % 
%
%  Set time span for integral for the ODE solver
dt_plot = 7*24*3600;
t_span = (0:dt_plot:6*dt_plot);
%
options_ode_z = odeset('RelTol',1e-8,'AbsTol',1e-9,'MaxStep',600);  
options_ode_r = odeset('RelTol',1e-8,'AbsTol',1e-9,'MaxStep',600);    
%
% Set up the time derivative and Jacobian function handles for the ODE solver
Temp_axial_t = @(t, T)Temp_axial_transient(t, T, hay, air);
Temp_radial_t = @(t, T)Temp_radial_transient(t, T, hay, air);
%
%  Run ODE integrator to get hay temperature variation with time
[hay.t, hay.T_z_t] = ode45(Temp_axial_t, t_span, T_z_t_0,options_ode_z);   % 
[hay.t, hay.T_r_t] = ode45(Temp_radial_t, t_span, T_r_t_0,options_ode_r);    %
%
%  Plot the hay temperature axial profile for each requested time increment 
figure;
plot(hay.z,hay.T_z_t(1,:), hay.z,hay.T_z_t(2,:), hay.z,hay.T_z_t(3,:), hay.z,hay.T_z_t(4,:), ...
       hay.z,hay.T_z_t(5,:), hay.z,hay.T_z_t(6,:), hay.z,hay.T_z_t(7,:),'LineWidth',2);
set(gca,'FontSize',14);
xlabel('width from center plance [m]','Fontsize',14);
ylabel('Temperature, [K]','Fontsize',14);
title('T(z,t) with T-Dependent Heat Generation')
legend(strcat('day ', num2str(hay.t(1)/(24*3600))), strcat('day ', num2str(hay.t(2)/(24*3600))), ...
           strcat('day ', num2str(hay.t(3)/(24*3600))), strcat('day ', num2str(hay.t(4)/(24*3600))), ...
           strcat('day ', num2str(hay.t(5)/(24*3600))), strcat('day ', num2str(hay.t(6)/(24*3600))), ...
           strcat('day ', num2str(hay.t(7)/(24*3600))));
%
%  Plot the hay temperature radial profile for each requested time increment
figure;
plot(hay.r,hay.T_r_t(1,:), hay.r,hay.T_r_t(2,:), hay.r,hay.T_r_t(3,:), hay.r,hay.T_r_t(4,:), ...
       hay.r,hay.T_r_t(5,:), hay.r,hay.T_r_t(6,:), hay.r,hay.T_r_t(7,:),'LineWidth',2); 
set(gca,'FontSize',14);
xlabel('Radius [m]','Fontsize',14);
ylabel('Temperature, [K]','Fontsize',14);
title('T(r,t) with T-Dependent Heat Generation')
legend(strcat('day ', num2str(hay.t(1)/(24*3600))), strcat('day ', num2str(hay.t(2)/(24*3600))), ...
           strcat('day ', num2str(hay.t(3)/(24*3600))), strcat('day ', num2str(hay.t(4)/(24*3600))), ...
           strcat('day ', num2str(hay.t(5)/(24*3600))), strcat('day ', num2str(hay.t(6)/(24*3600))), ...
           strcat('day ', num2str(hay.t(7)/(24*3600))));
end 
%
%% Function for calculating time derivatives for axial temperature profile
function [dT_dt]  = Temp_axial_transient(t, T, hay, air)
% Divide Residual by density and specific heat
%
[res] = Temp_axial(T, hay, air);
dT_dt = res/(hay.c_v*hay.rho_b);   % Derivates of temperature at each axial node [K s^-1]
end
%
%% Function for calculating time derivatives for radial temperature profile
function [dT_dt]  = Temp_radial_transient(t, T, hay, air)
%  Divide Residual by density and specific heat
%
[res] = Temp_radial(T, hay, air);
dT_dt = res/(hay.c_v*hay.rho_b);    % Derivates of temperature at each radial node [K s^-1]
end
%
%% Function for calculating residuals of axial temperature profile
function [res]  = Temp_axial(T, hay, air)
%
%  Write comments showing numerical discretization of heat conduction equation for interior points
%   Start with HDE rho*cv*dT/dt = qin"*A + qout"*A qdot*V
%   Use forward half step finite volume expression for heat leaving
%   Use backward half step finite volume expression for heat entering
%   Circular areas do not change with Z, they cancel w/ V leaving 1/dz

%  Write comments showing discretization of z = 0  boundary 
%   At the adiabatic condition, Areas and volumes cancel out to 1/dz
%   There is no qin" as there is no backward half step
%   Sub in Forward half step expression for qout" yields next step minus
%   current step terms

%  Write comments showing discretization of z = z_max  boundary
%  Heat enters and leaves through through identical areas, volume is over 2
%  because there is no forward step as we are at zmax.
%  Sub in axial step expression for qin" and Thermal circuit for qout" to
%  get a uniqe boundary condition

%  Initialize energy balance residuals for all of the axial nodes.
hay.q_reac_z = @(T) hay.A_reac * exp(-hay.T_a_reac / T);
res = -999*ones(length(T),1);
%
% Set up residual for z = 0 boundary cell 
res(1) = 2*hay.k_b*(T(2)-T(1))/hay.dz^2 + hay.q_reac_z(T(1)); 

%
% Set up residuals for interior cells 
for i_z = 2:hay.n_z - 1
    res(i_z) = -hay.k_b*(T(i_z)-T(i_z-1))/(hay.dz^2) + hay.k_b*(T(i_z+1)- T(i_z))/(hay.dz^2) +hay.q_reac_z(T(i_z));  % ???  
end
% Set up residual for z = 0.5*z_max boundary cell 
res(end) = -2*hay.k_b*(T(end)-T(end-1))/(hay.dz^2) + 2*(air.T - T(end)) / (hay.R_t_ext*hay.dz) + hay.q_reac_z(T(end));  % ???  
end
%
%% Function for calculating residuals of radial temperature profile
function [res]  = Temp_radial(T, hay, air)
%
%  Write comments showing numerical discretization of heat conduction equation for interior points  
%  %Start with HDE rho*cv*dT/dt = qin"*A + qout"*A qdot*V
%   Use forward half step finite volume expression for heat leaving
%   Use backward half step finite volume expression for heat entering
%   Areas and volumes areas change with r, so use backward (r-dr) and forward area
%   expressions (r+dr) for q"in and qout" respectively
%   Include 0.5 to account for half step area

%  Write comments showing discretization of r = 0  boundary  
%  Adiabatic Centerline so there is no qin"
%  susbstitute dr/2, dphi/2,  for volume as we are looking at a half
%  cylinder, and dr/2 for area in our half step forward
% Substitute q"out for forward half step finite volume expression
%
%  Write comments showing discretization of r = r_max  boundary
%  Heat enters and leaves through through hay-plastic interface
%  Substitute in backward half step expression for qout" and Backward area
%  half step for Ain. 
%  Sub in Thermal circuit for qout" and max area
%  Volume is confined between the max and a half step back
%

%  Initialize energy balance residuals for all of the radial nodes.
hay.q_reac_z = @(T) hay.A_reac * exp( -hay.T_a_reac / T);
res = -999*ones(length(T),1);
%
% Set up residual for r = 0 boundary 
% res(1) = 8*hay.k_b*(T(2)-T(1))/(hay.dr^3) + hay.q_reac_z(T(1));   % ???
res(1) = 4*hay.k_b*(T(2)-T(1))/(hay.dr^2) + hay.q_reac_z(T(1));
%
% Set up residuals for interior points 
for i_r = 2:hay.n_r - 1
    res(i_r) = -0.5*hay.k_b/hay.r(i_r)*(hay.r(i_r-1)+hay.r(i_r))*(T(i_r)-T(i_r-1))/(hay.dr^2) + hay.k_b/hay.r(i_r)*(0.5*(hay.r(i_r)+hay.r(i_r+1)))*(T(i_r+1)-T(i_r))/(hay.dr^2) + hay.q_reac_z(T(i_r));
end
%
% Set up residual for r = r_max boundary 
res(end)  = (-hay.k_b*(hay.r_max - 0.5*hay.dr)*2*(T(end)-T(end-1))/hay.dr + 2*hay.r_max*(air.T-T(end))/hay.R_t_ext)/(hay.r_max*hay.dr - 0.25*hay.dr^2) + hay.q_reac_z(T(end));
end