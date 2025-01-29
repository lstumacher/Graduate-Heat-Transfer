function[pipe, fluid, air, mesh]  = MEGN571_HW5_F2024_starter()
%
%%  MEGN571 Module 5 Assignment Fall 2024 starter code:
%       by Prof. Greg Jackson, Colorado School of Mines
%
%  See Module 5 assignment write-up for the problem description
%  Solution to thermally developing laminar pipe flow
%  
%% Initialize properties of the pipe, fluid inlet, and the external air flows.
close all;
close all;
%
const.g = 9.81;   % Earth's gravity acceleration [m/s^2]
%
pipe.L_z = 4.0;	% length of pipe [m]
pipe.r_int = 0.04;	% internal radius of pipe [m]
pipe.r_ext = 0.05;	% external radius of pipe [m]
pipe.A_c = pi*pipe.r_int^2;     % cross-sectional area for flow [m^2]
%
fluid.T_init = 180+273.15;	% inlet temperature of water [K]
fluid.P_init = 2e6;   % inlet pressure of water [Pa]
fluid.mdot = 0.02;	% mass flow through pipe [kg/s]
%
%    ???
fluid.c_p = py.CoolProp.CoolProp.PropsSI('C', 'T', fluid.T_init, 'P', fluid.P_init, 'Water');; 	% specific heat capacity of fluid [J/kg-K]
fluid.rho = py.CoolProp.CoolProp.PropsSI('D', 'T', fluid.T_init, 'P', fluid.P_init, 'Water');; 	% density of fluid [kg/m^3]
fluid.k = py.CoolProp.CoolProp.PropsSI('L', 'T', fluid.T_init, 'P', fluid.P_init, 'Water');;     % thermal conductivity of fluid [W/m-K]
fluid.mu = py.CoolProp.CoolProp.PropsSI('V', 'T', fluid.T_init, 'P', fluid.P_init, 'Water'); % fluid viscosity [kg/m-s]
fluid.u_m = fluid.mdot/(fluid.rho*pipe.A_c);  % mean velocity through pipe [m/s]
%
air.obj = Solution('air.yaml','air');  % Cantera object to get air properties 
air.P = 100000;                % presure of surrounding air [Pa]
air.T = 10+273.15;	                    % temperature of surrounding air [K]
%
pipe.T_s_a = 0.5*(fluid.T_init+air.T);   % constant pipe wall temperature for parts a and b [deg. C]
pipe.k = 16;    % thermal conductivity of pipe [W/m-K]
%
% Calculate relevant dimensionless properties for the internal fluid flow
fluid.Pr = fluid.mu*fluid.c_p/(fluid.k);                          % fluid Prandtl number [--]
fluid.Re_D = fluid.rho*fluid.u_m*pipe.r_int*2/fluid.mu;  % Reynolds number of internal flow [--]
fluid.Gz_inv = pipe.L_z/(fluid.Pr*fluid.Re_D);   % inverse of Graetz number [--]
%  fluid.f =  ???          % friction factor (Petukhov) [--]
%
%% Set up numerical mesh, radial slice cross-sectional areas, and fluid axial velocity radial profile
%
%% Set up numerical mesh, radial slice cross-sectional areas, and fluid axial velocity radial profile
mesh.dr = pipe.r_int/100;               % mesh radial discertization [m]
mesh.dz = pipe.L_z/100;                % mesh axial discertization [m]
mesh.r = (0:mesh.dr:pipe.r_int);      % mesh radial locations [m]
mesh.z = (0.0:mesh.dz:pipe.L_z)';    % mesh axial locations [m]
mesh.n_r = length(mesh.r);               % number of radial nodes [--]
mesh.n_z = length(mesh.z);             % number of axial nodes [--]
%
mesh.A_r = zeros(1,mesh.n_r);    % areas for radial conduction on the outer boundary of cell i_r [m^2]
mesh.A_z = zeros(1,mesh.n_r);   % cross-sectional area for the z-direction flow at cell i_r [m^2]
fluid.u_z = zeros(1,mesh.n_r);      % axial velocity radial profile for fully developed flow.
%
mesh.A_z(1,1) =  pi*((mesh.r(1)+mesh.dr/2)^2);  
mesh.A_r(1,1) =  mesh.dz*2*pi*(mesh.r(1)+mesh.dr/2);  
fluid.u_z(1,1) =  2*fluid.u_m;
for i_r = 2:mesh.n_r-1
  mesh.A_z(1,i_r) =  pi*((mesh.r(i_r)+mesh.dr/2).^2-(mesh.r(i_r)-mesh.dr/2).^2);  
  mesh.A_r(1,i_r) =  mesh.dz*2*pi*(mesh.r(i_r)+mesh.dr/2); 
  fluid.u_z(1,i_r) =  2*fluid.u_m*(1-(mesh.r(i_r)/pipe.r_int)^2); 
end
 mesh.A_z(1,mesh.n_r) =  pi*((mesh.r(mesh.n_r))^2-(mesh.r(mesh.n_r)-mesh.dr/2)^2); 
 mesh.A_r(1,mesh.n_r) =  mesh.dz*2*pi*(pipe.r_int); 
 fluid.u_z(1,mesh.n_r) =  fluid.u_m;%0; 
%
%% Part a) Find Nusselt number correlation for thermally developing flow
%         as a function of z using Graetz number as a function 
fluid.Gz_D = fluid.Re_D*fluid.Pr*2*pipe.r_int./mesh.z;
fluid.Nu_D_corr = 3.657+0.2362*(fluid.Gz_D.^0.488).*exp(-57.2./fluid.Gz_D);
%
% Make plot of correlation Nusselt number with respect to z 
% from Graetz number correlation
figure(1)
semilogx(mesh.z, fluid.Nu_D_corr,'LineWidth',2);
axis([0 pipe.L_z 0 50])
xlabel('axial distance, z [m]','FontSize',16)
ylabel('Correlation Nusselt Number -- Nu_z','FontSize',16)
%
% Check the mean velocity to see if it matches fluid.u_m
u_m = sum(fluid.u_z.*mesh.A_z)/pipe.A_c;
%
%% Part b)  Solve for radial temperature profiles for all z with constant wall T
options = odeset('RelTol',1e-6);
%
%  Initialize T at z = 0
T_0 = fluid.T_init*ones(mesh.n_r,1);
T_0(end) = pipe.T_s_a;
is_T_s_const = true;
%
% Integrate dTdz function with respect to z
z_span = mesh.z;
[z, T] = ode45(@(z,T) dTdz(z, T, mesh, pipe, fluid, air, const, is_T_s_const), z_span, T_0, options);
fluid.T_b = T;
%
% Set the indices for the z_plot and plot the radial profile at all z_plot value 
z_plot = pipe.L_z*[0; 0.2; 0.4; 0.6; 0.8; 1.0];
i_z_plot = zeros(length(z_plot),1);
for i_z = 1:length(z_plot)
    i_z_plot(i_z,1) = find(mesh.z == z_plot(i_z));
end
%
figure(2);
plot(mesh.r, fluid.T_b(i_z_plot,:)-273.15,'LineWidth',2);
axis([0 pipe.r_int 90 200])
xlabel('radius, r [m]','FontSize',16)
ylabel('Temperatures [°C]','FontSize',16)
legend(num2str(z_plot),'Location','best')
%
%  Set up variables for finding mean fluid temperature and Nusselt number, and
%  plot them both as s a function of z
fluid.T_m_b = zeros(mesh.n_z,1);
fluid.Nu_z_num = zeros(mesh.n_z,1);
for i_z =1:mesh.n_z
       integrand = 2*pi*mesh.r.*fluid.u_z.*fluid.T_b(i_z,:);
       fluid.T_m_b(i_z,1) = trapz(mesh.r, integrand)./(u_m.*pi.*(pipe.r_int.^2));
       qflux = -fluid.k*(fluid.T_b(i_z,mesh.n_r-1)-pipe.T_s_a)/mesh.dr;
       fluid.Nu_z_num_b(i_z,1) = qflux*2*pipe.r_int/(fluid.k*(pipe.T_s_a-fluid.T_m_b(i_z,1)));
       
end

%
figure(3);
plot(mesh.z, fluid.T_m_b-273.15, mesh.z, fluid.T_b(:,mesh.n_r)-273.15,'LineWidth',2);
axis([0 pipe.L_z 0 200])
xlabel('axial distance, z [m]','FontSize',16)
ylabel('Mean and Wall Temperatures [°C]','FontSize',16)
legend('Mean fluid','Wall','Location','best')
%
figure(4);
semilogx(mesh.z, fluid.Nu_z_num_b, mesh.z, fluid.Nu_D_corr,'LineWidth',2);
axis([0 pipe.L_z 0 50])
xlabel('axial distance, z [m]','FontSize',16)
ylabel('Nusselt Number -- Nu_z','FontSize',16)
legend('Numerical','Correlation')
%
%% Part c)  Solve for radial temperature profile for all z with external thermal resistance to air
%
 % Initialize T at z = 0
T_0 = fluid.T_init*ones(mesh.n_r,1);
T_0(end) = pipe.T_s_a;
is_T_s_const = false;
z_span = mesh.z;
[z, T] = ode45(@(z,T) dTdz(z, T, mesh, pipe, fluid, air, const, is_T_s_const), z_span, T_0, options);
fluid.T_c = T;
%
% Plot the radial profile at all z_plot value 
figure(5);
plot(mesh.r, T(i_z_plot,:)-273.15,'LineWidth',2);
axis([0 pipe.r_int 0 200])
xlabel('radius, r [m]','FontSize',16)
ylabel('Temperatures [°C]','FontSize',16)
legend(num2str(z_plot),'Location','best')
%
%  Set up variables for finding mean fluid temperature and Nusselt number, and
%  plot them both as s a function of z
fluid.T_m_c = zeros(mesh.n_z,1);
fluid.Nu_z_num_c = zeros(mesh.n_z,1);
for i_z =1:mesh.n_z

       integrand = 2*pi*mesh.r.*fluid.u_z.*fluid.T_c(i_z,:);
       fluid.T_m_c(i_z,1) = trapz(mesh.r, integrand)./(u_m.*pi.*(pipe.r_int.^2));
       qflux = -fluid.k*(fluid.T_c(i_z,mesh.n_r-1)-pipe.T_s_a)/mesh.dr;
       fluid.Nu_z_num_c(i_z,1) = qflux*2*pipe.r_int/(fluid.k*(pipe.T_s_a-fluid.T_m_c(i_z,1)));
end
%
figure(6);
plot(mesh.z, fluid.T_m_c-273.15, mesh.z, fluid.T_c(:,mesh.n_r)-273.15,'LineWidth',2);
axis([0 pipe.L_z 0 200])
xlabel('axial distance, z [m]','FontSize',16)
ylabel('Mean and Wall Temperatures [°C]','FontSize',16)
legend('Mean fluid','Wall','Location','best')
%
figure(7);
semilogx(mesh.z, fluid.Nu_z_num_c, mesh.z, fluid.Nu_D_corr,'LineWidth',2);
axis([0 pipe.L_z 0 50])
xlabel('axial distance, z [m]','FontSize',16)
ylabel('Nusselt Number -- Nu_z','FontSize',16)
legend('Numerical','Correlation')
 end
% %
%%  Function that calculates dTdz for an ode integrator
function [dTdz] = dTdz(z, T, mesh, pipe, fluid, air, const, is_T_s_const)
% Inputs
% z - position (m)
% T - vector of temperatures (K)
% mesh, pipe, fluid, air are structures with properties for solving problem
%  is_T_s_const
%
% Calculate dTdz - rate of temperature change in z for each node (K/m)
dTdz = zeros(mesh.n_r,1);
dTdz(1) =  -mesh.A_r(1)*fluid.k*((T(1)-T(2))/mesh.dr)/(fluid.rho*fluid.u_z(1)*fluid.c_p*mesh.A_z(1)*mesh.dz);
for i_r = 2:(mesh.n_r-1)
    dTdz(i_r) = (mesh.A_r(i_r-1)*fluid.k*((T(i_r-1)-T(i_r))/mesh.dr)...
       -(mesh.A_r(i_r)*fluid.k)*((T(i_r)-T(i_r+1))/mesh.dr))...
       /(fluid.rho*fluid.u_z(i_r)*fluid.c_p*mesh.A_z(i_r)*mesh.dz) ;
end
%
%  Calculate different wall cell derivative based on whether T is constant.
if is_T_s_const
   dTdz(mesh.n_r) = 0; %% 2*pi*(pipe.r_int-mesh.dr/2)*fluid.k*((T(i_z,mesh.n_r-1)-T(i_z,mesh.n_r))/mesh.dr); 
else
    %
    %  Find external heat transfer coefficient for external thermal
    %  resistance using Churchill and Chu eqn. for natural convection Nu_D
    T_film = 0.5*(air.T + T(mesh.n_r));
    %disp(['T_film = ', num2str(T_film)]);
    set(air.obj, 'P', air.P, 'T', T_film)
    air.rho = density(air.obj);
    air.mu = viscosity(air.obj);
    air.k = thermalConductivity(air.obj);
    air.cp = cp_mass(air.obj);
    air.Pr = air.mu*air.cp/air.k;
    %disp(air.Pr)
    air.Beta = (1/T_film);
    Ra = const.g *  (1/T_film) * (T(mesh.n_r)-air.T) * ((2*pipe.r_ext)^3) * (air.rho^2)*(air.cp) / (air.mu * air.k);
    %RaP = const.g*air.Beta*(T(mesh.n_r)-air.T)*(pipe.r_ext*2)^3*air.rho^2*air.cp/(air.mu*air.k);
    NuD = 0.36 + 0.518*(Ra^0.25)/((1 + (0.559/air.Pr)^(9/16))^(4/9));
    h_T_air = NuD*air.k/(2*pipe.r_ext);
    R_th_ext = 1/(h_T_air) + log(pipe.r_ext/pipe.r_int)/(pipe.k);
    % -k dT/dz = h(T-T)
    %  %  Set up derivative for boundary cell ignoring pipe axial conduction
    dTdz(mesh.n_r) = ((mesh.A_r(mesh.n_r - 1) * fluid.k / mesh.dz) * (T(mesh.n_r - 1) - T(mesh.n_r)) / mesh.dr ...
                         +(mesh.A_r(mesh.n_r) / mesh.dz) * (air.T - T(mesh.n_r)) / R_th_ext) ...
                         / (mesh.A_z(mesh.n_r) * fluid.rho * fluid.u_z(mesh.n_r) * fluid.c_p);

end 
end
