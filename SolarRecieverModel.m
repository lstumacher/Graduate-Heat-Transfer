function[pipe, salt, air, mesh]  = MEGN571_ProjectIteration3()
%
%%  MEGN571 Module 5 Assignment Fall 2024 starter code:
%       by Prof. Greg Jackson, Colorado School of Mines
%
%  See Module 5 assignment write-up for the problem description
%  Solution to thermally developing laminar pipe flow
%  
%% Initialize properties of the pipe, salt inlet, and the external air flows.
close all;
close all;
%
const.g = 9.81;% Earth's gravity acceleration [m/s^2]
const.sigma = 5.67*10^-8; 
%
pipe.L_z = 30;	% length of pipe [m]
pipe.r_int = 7.89835;	% internal radius of pipe [m]
pipe.r_ad = 7.87335; % radius of internal pipe
pipe.r_ext = 7.9;	% external radius of pipe [m]
pipe.A_c = pi*(pipe.r_int^2-pipe.r_ad^2);     % cross-sectional area for flow [m^2]
pipe.rho = 9000;
pipe.cp = 200;
%
salt.T_init = 290+273.15;	% inlet temperature of water [K]
salt.T_outlet = 565+273.15;	% Outlet temperature of water [K]
salt.T_av = (salt.T_init + salt.T_outlet)/2;

%salt.P_init = 2e6;   % inlet pressure of water [Pa]
salt.mdot = 67;	% mass flow through pipe [kg/s]
salt.mdot_values = [40, 67, 85]; % Mass flow rates for parametric study [kg/s]
%
%    ???
salt.c_p = 1442.3; 	% specific heat capacity of salt [J/kg-K]
salt.rho = 2085.5; 	% density of salt [kg/m^3]
salt.k = 0.4407;     % thermal conductivity of salt [W/m-K]
salt.mu = 1.529*(10^(-3)); % salt viscosity [kg/m-s]
salt.u_m = salt.mdot/(salt.rho*pipe.A_c);  % mean velocity through pipe [m/s]
%
air.P = 100000;                % presure of surrounding air [Pa]
air.T = 20+273.15;	                    % temperature of surrounding air [K]
air.h = 32; %Convective heat transfer coefficient
air.q_flux = 445000;
air.q_flux_t = sort([431835.4863 608998.7627 694416.7709 742662.1275 770343.8894 770343.8894 777462.0568 771925.7044 760062.0921 738707.59 697580.4009 627980.5423 506180.7897 244390.4126]); %Solar Flux
%air.q_flux_t = [136800, 196800];
%
pipe.T_s_a = 565+273.15;   % constant pipe wall temperature for parts a and b [deg. C]
pipe.k = 400;% thermal conductivity of pipe [W/m-K]
pipe.alpha = 0.95; %absorptivity
pipe.eps = 0.89; %emissivity

%
% Calculate relevant dimensionless properties for the internal salt flow
salt.Pr = salt.mu*salt.c_p/salt.k;                       % salt Prandtl number [--]
salt.Re_D = salt.rho*salt.u_m*(pipe.r_int*2-pipe.r_ad*2)/salt.mu;  % Reynolds number of internal flow [--]
%salt.Gz_inv = =  ???   % inverse of Graetz number [--]
%salt.f =  ???          % friction factor (Petukhov) [--]
%
 
%% Set up numerical mesh, radial slice cross-sectional areas, and salt axial velocity radial profile
mesh.dr = (pipe.r_int-pipe.r_ad)/100;               % mesh radial discertization [m]
mesh.dz = pipe.L_z/100;                % mesh axial discertization [m]
mesh.r = (pipe.r_ad:mesh.dr:pipe.r_int);      % mesh radial locations [m]
mesh.r = [mesh.r, pipe.r_ext];
mesh.z = (0.0:mesh.dz:pipe.L_z)';   % mesh axial locations [m]
mesh.n_r = length(mesh.r);               % number of radial nodes [--]
mesh.n_z = length(mesh.z);             % number of axial nodes [--]
%
mesh.A_r = zeros(1,mesh.n_r);    % areas for radial conduction on the outer boundary of cell i_r [m^2]
mesh.A_z = zeros(1,mesh.n_r);   % cross-sectional area for the z-direction flow at cell i_r [m^2]
salt.u_z = zeros(1,mesh.n_r-1);      % axial velocity radial profile for fully developed flow.
%
mesh.A_z(1,1) =  pi*((mesh.r(1)+mesh.dr/2)^2 - (mesh.r(1)-mesh.dr/2)^2) ;  
mesh.A_r(1,1) =  2*pi*mesh.dz*(mesh.r(1)+mesh.dr/2);  
salt.u_z(1,1) =  0;%2*salt.u_m;
for i_r = 2:mesh.n_r-2

  mesh.A_z(1,i_r) =  pi*((mesh.r(i_r)+mesh.dr/2)^2-(mesh.r(i_r)-mesh.dr/2)^2);  
  mesh.A_r(1,i_r) =  2*pi*mesh.dz*(mesh.r(i_r)+mesh.dr/2); 
  salt.u_z(1,i_r) =  1.5*salt.u_m * (1 - ((2*mesh.r(i_r)-(pipe.r_ad+pipe.r_int))/(pipe.r_int-pipe.r_ad)).^2); %2 * salt.u_m * (1 - ((mesh.r(i_r) - pipe.r_ad) * (mesh.r(i_r) - pipe.r_int)) / ((pipe.r_ad - (pipe.r_ad + pipe.r_int) / 2) * (pipe.r_int - (pipe.r_ad + pipe.r_int) / 2)));%2*salt.u_m*(1-(mesh.r(i_r)/pipe.r_int)^2);
end
mesh.A_z(1,mesh.n_r-1) =  pi*((mesh.r(mesh.n_r-1)+mesh.dr/2)^2 - (mesh.r(mesh.n_r-1)-mesh.dr/2)^2) ;   
mesh.A_r(1,mesh.n_r-1) =  2*pi*mesh.dz*(mesh.r(mesh.n_r-1)+mesh.dr/2); 
salt.u_z(1,mesh.n_r-1) =  0;

%Adding another node for the outer surface 
mesh.A_z(1,mesh.n_r) =  pi*((mesh.r(end)+mesh.dr/2)^2 - (mesh.r(end)-mesh.dr/2)^2) ;   
mesh.A_r(1,mesh.n_r) =  2*pi*mesh.dz*(mesh.r(end)+mesh.dr/2); 
% salt.u_z(1,mesh.n_r) =  0; %2*salt.u_m*(1-(mesh.r(mesh.n_r-1)/pipe.r_int)^2);
% mass_flow_rate = salt.rho * salt.u_z .* (2 * pi * mesh.r);
% mass_flow_rate_value = trapz(mesh.r, mass_flow_rate);
%
%% Part a) Find Nusselt number correlation for thermally developing flow
%         as a function of z using Graetz number as a function 
salt.Gz_D = salt.Re_D*salt.Pr*2*(pipe.r_int-pipe.r_ad)/pipe.L_z;
salt.Nu_D_corr = 3.66 + 0.2362*(salt.Gz_D^(0.488))*exp(-57.2/salt.Gz_D);
salt.h = salt.Nu_D_corr*salt.k/(2*(pipe.r_int-pipe.r_ad));
%
% Check the mean velocity to see if it matches salt.u_m
%u_m = sum(salt.u_z.*mesh.A_z(1,1:mesh.n_r-1))/pipe.A_c;
%
%%  Solve for radial temperature profiles for all z 
options = odeset('RelTol',1e-6);
%
%  Initialize T at z = 0
T_0 = salt.T_init*ones(mesh.n_r,1);
%T_0(end) = (salt.T_init + air.T)/2;
T_wall_0 = T_0(end);
%
% Integrate dTdz function with respect to z
z_span = mesh.z;
[z, T] = ode15s(@(z,T) dTdz( z, T, mesh, pipe, salt, air, const),z_span,T_0,options);
salt.T_b = T;

% Set the indices for the z_plot and plot the radial profile at all z_plot value 
z_plot = pipe.L_z*[0; 0.2; 0.4; 0.6; 0.8; 1.0];
i_z_plot = zeros(length(z_plot),1);
for i_z = 1:length(z_plot)
    i_z_plot(i_z,1) = find(mesh.z == z_plot(i_z));
end


figure;
plot(mesh.r, salt.T_b(i_z_plot,:)-273.15,'LineWidth',2);
%axis([pipe.r_ad pipe.r_int])
xlabel('radius, r [m]','FontSize',16)
ylabel('Temperatures [°C]','FontSize',16)
legend(num2str(z_plot),'Location','best')
%
for i_z = 1:mesh.n_z
    salt.T_m(i_z,1) = sum((mesh.A_z(1,1:end-1).*salt.u_z).* salt.T_b(i_z,1:end-1))/(pipe.A_c*salt.u_m);
end

figure(2)
plot(mesh.z, salt.T_m -273.15, mesh.z,salt.T_b(:,mesh.n_r-1)-273.15,'LineWidth',2)

xlabel('axial distance, z [m]','FontSize',16)
ylabel('Mean and Wall Temperature [C]','FontSize',16)
legend('Mean fluid','Wall','Location','best')

disp(salt.T_m(1,1)-273.15)
disp(salt.T_m(end,1)-273.15)

%% Parametric study
air.q_solar_var = [450000 608998.7627 694416.7709 742662.1275 770343.8894 770343.8894 777462.0568 771925.7044 760062.0921 738707.59 697580.4009 627980.5423 506180.7897 244390.4126];
air.q_solar_var = sort(air.q_solar_var);
salt.T_m_flux = zeros(mesh.n_z, length(air.q_solar_var));
salt.mdot_var = [40 67 85];
for i = 1:length(salt.mdot_var)
salt.mdot = salt.mdot_var(i);
salt.u_m = salt.mdot/(salt.rho*pipe.A_c);
for i_r = 2:mesh.n_r-2
salt.u_z(1,i_r) = 1.5*salt.u_m * (1 - ((2*mesh.r(i_r)-(pipe.r_ad+pipe.r_int))/(pipe.r_int-pipe.r_ad)).^2); %2 * salt.u_m * (1 - ((mesh.r(i_r) - pipe.r_ad) * (mesh.r(i_r) - pipe.r_int)) / ((pipe.r_ad - (pipe.r_ad + pipe.r_int) / 2) * (pipe.r_int - (pipe.r_ad + pipe.r_int) / 2)));%2*salt.u_m*(1-(mesh.r(i_r)/pipe.r_int)^2);
end
for i_flux = 1:length(air.q_solar_var)
air.q_flux = air.q_solar_var(i_flux);
[~, T_flux] = ode15s(@(z,T) dTdz(z, T, mesh, pipe, salt, air, const), z_span, T_0, options);
for i_z = 1:mesh.n_z
salt.T_m_flux(i_z, i_flux) = sum((mesh.A_z(1,1:end-1).*salt.u_z).* T_flux(i_z,1:end-1))/(pipe.A_c*salt.u_m);
end
end
figure(3);
hold on
plot(air.q_solar_var, salt.T_m_flux(end-1,:) - 273.15, '-o');
xlabel('Solar Flux, q [W/m^2]','FontSize',16);
ylabel('Mean Fluid Outlet Temperature [\circC]','FontSize',16);
legend({'$\dot{m}$ = 40 kg/s','$\dot{m}$ = 67 kg/s','$\dot{m}$ = 85 kg/s'}, 'Interpreter', 'latex', 'Location', 'best');
title('Parametric Study: Effect of Solar Flux on Outlet Mean Fluid Temperature');
end

%%
for i_mdot = 1:length(salt.mdot_values)
        % Update salt and air properties for this iteration
        salt.mdot = salt.mdot_values(i_mdot);
        salt.u_m = salt.mdot/(salt.rho*pipe.A_c); % mean velocity through pipe [m/s]
    for i_r = 2:mesh.n_r-2
        salt.u_z(1,i_r) = 1.5*salt.u_m * (1 - ((2*mesh.r(i_r)-(pipe.r_ad+pipe.r_int))/(pipe.r_int-pipe.r_ad)).^2); %2 * salt.u_m * (1 - ((mesh.r(i_r) - pipe.r_ad) * (mesh.r(i_r) - pipe.r_int)) / ((pipe.r_ad - (pipe.r_ad + pipe.r_int) / 2) * (pipe.r_int - (pipe.r_ad + pipe.r_int) / 2)));%2*salt.u_m*(1-(mesh.r(i_r)/pipe.r_int)^2);
    end

    for i_qflux = 1:length(air.q_solar_var)
        
        air.q_flux = air.q_solar_var(i_qflux);
        [z, T] = ode15s(@(z,T) dTdz(z, T, mesh, pipe, salt, air, const), z_span, T_0, options);
        for i_z = 1:mesh.n_z
            salt.T_m(i_z, i_qflux) = sum((mesh.A_z(1,1:end-1).*salt.u_z).* T(i_z,1:end-1))/(pipe.A_c*salt.u_m);
        end

                % Calculate final mean temperature
        % T_final = T(end, :);
        % salt.T_m_final = sum(T_final(1:end-1).*mesh.r(1:end-1)) / sum(mesh.r(1:end-1));
        % final_temperatures(i_mdot, i_qflux) = salt.T_m_final - 273.15;
        % Calculate mean temperature at each z
        % salt.T_m = zeros(mesh.n_z, 1);
        % for i_z = 1:mesh.n_z
        %     salt.T_m(i_z,1) = sum(mesh.A_z(1:end-1).*salt.u_z.*T(i_z,1:end-1)) / (pipe.A_c*salt.u_m);
        %     %salt.T_m(i_z, i_qflux) = sum((mesh.A_z(1,1:end-1).*salt.u_z).* T(i_z,1:end-1))/(pipe.A_c*salt.u_m);
        % end

        % T_final = T(end, :);
        % salt.T_m_final = sum(T_final(1:end-1).*mesh.r(1:end-1)) / sum(mesh.r(1:end-1));
        % final_temperatures(i_mdot, i_qflux) = salt.T_m_final - 273.15;
    end
figure(4);
hold on
plot(air.q_solar_var, salt.T_m(end-1,:) - 273.15, '-o');
xlabel('Solar Flux, q [W/m^2]','FontSize',16);
ylabel('Mean Fluid Outlet Temperature [\circC]','FontSize',16);
legend({'$\dot{m}$ = 40 kg/s','$\dot{m}$ = 67 kg/s','$\dot{m}$ = 85 kg/s'}, 'Interpreter', 'latex', 'Location', 'best');
title('Parametric Study: Effect of Solar Flux on Outlet Mean Fluid Temperature');
end
% Format the plot
% 
% figure;
% hold on;
% for i_mdot = 1:length(salt.mdot_values)
%     plot(air.q_flux_t, final_temperatures(i_mdot, :), '-o', 'DisplayName', sprintf('mdot = %.1f kg/s', salt.mdot_values(i_mdot)));
% end
% xlabel('Solar Flux [W/m^2]', 'FontSize', 16);
% ylabel('Final Mean Temperature [\circC]', 'FontSize', 16);
% title('Final Mean Temperature vs. Heat Flux', 'FontSize', 16);
% legend('show', 'Location', 'best');
% grid on;


%% Tube Size Parametric Study with mass

%Standard
air.q_flux = 445000;
salt.mdot_var = [40 67 85];

%Variable
pipe.L_z_var = 5:1:50;	% length of pipe [m]
pipe.r_int = 7.89835;	% internal radius of pipe [m]
pipe.r_ad = 7.87335; % radius of internal pipe
pipe.r_ext_var = [7.899, 7.9, 8];	% external radius of pipe [m]
for i_mdot = 1:length(salt.mdot_values)
        % Update salt and air properties for this iteration
        salt.mdot = salt.mdot_values(i_mdot);
        salt.u_m = salt.mdot/(salt.rho*pipe.A_c); % mean velocity through pipe [m/s]
for i_L_z = 1:length(pipe.L_z_var)
    pipe.L_z = pipe.L_z_var(i_L_z);
    pipe.A_c = pi*(pipe.r_int^2-pipe.r_ad^2);     % cross-sectional area for flow [m^2]
    pipe.rho = 9000;
    pipe.cp = 200;
    mesh.dr = (pipe.r_int-pipe.r_ad)/100;               % mesh radial discertization [m]
    mesh.dz = pipe.L_z/100;                % mesh axial discertization [m]
    mesh.r = (pipe.r_ad:mesh.dr:pipe.r_int);      % mesh radial locations [m]
    mesh.r = [mesh.r, pipe.r_ext];
    mesh.z = (0.0:mesh.dz:pipe.L_z)';   % mesh axial locations [m]
    mesh.n_r = length(mesh.r);               % number of radial nodes [--]
    mesh.n_z = length(mesh.z);             % number of axial nodes [--]
    %
    mesh.A_r = zeros(1,mesh.n_r);    % areas for radial conduction on the outer boundary of cell i_r [m^2]
    mesh.A_z = zeros(1,mesh.n_r);   % cross-sectional area for the z-direction flow at cell i_r [m^2]
    salt.u_z = zeros(1,mesh.n_r-1);      % axial velocity radial profile for fully developed flow.
    %
    mesh.A_z(1,1) =  pi*((mesh.r(1)+mesh.dr/2)^2 - (mesh.r(1)-mesh.dr/2)^2) ;  
    mesh.A_r(1,1) =  2*pi*mesh.dz*(mesh.r(1)+mesh.dr/2);  
    salt.u_z(1,1) =  0;%2*salt.u_m;
    for i_r = 2:mesh.n_r-2
    
      mesh.A_z(1,i_r) =  pi*((mesh.r(i_r)+mesh.dr/2)^2-(mesh.r(i_r)-mesh.dr/2)^2);  
      mesh.A_r(1,i_r) =  2*pi*mesh.dz*(mesh.r(i_r)+mesh.dr/2); 
      salt.u_z(1,i_r) =  1.5*salt.u_m * (1 - ((2*mesh.r(i_r)-(pipe.r_ad+pipe.r_int))/(pipe.r_int-pipe.r_ad)).^2); %2 * salt.u_m * (1 - ((mesh.r(i_r) - pipe.r_ad) * (mesh.r(i_r) - pipe.r_int)) / ((pipe.r_ad - (pipe.r_ad + pipe.r_int) / 2) * (pipe.r_int - (pipe.r_ad + pipe.r_int) / 2)));%2*salt.u_m*(1-(mesh.r(i_r)/pipe.r_int)^2);
    end
    mesh.A_z(1,mesh.n_r-1) =  pi*((mesh.r(mesh.n_r-1)+mesh.dr/2)^2 - (mesh.r(mesh.n_r-1)-mesh.dr/2)^2) ;   
    mesh.A_r(1,mesh.n_r-1) =  2*pi*mesh.dz*(mesh.r(mesh.n_r-1)+mesh.dr/2); 
    salt.u_z(1,mesh.n_r-1) =  0;
    
    %Adding another node for the outer surface 
    mesh.A_z(1,mesh.n_r) =  pi*((mesh.r(end)+mesh.dr/2)^2 - (mesh.r(end)-mesh.dr/2)^2) ;   
    mesh.A_r(1,mesh.n_r) =  2*pi*mesh.dz*(mesh.r(end)+mesh.dr/2); 

        salt.u_m = salt.mdot/(salt.rho*pipe.A_c); % mean velocity through pipe [m/s]
        for i_r = 2:mesh.n_r-2
                salt.u_z(1,i_r) =  1.5*salt.u_m * (1 - ((2*mesh.r(i_r)-(pipe.r_ad+pipe.r_int))/(pipe.r_int-pipe.r_ad)).^2); %2 * salt.u_m * (1 - ((mesh.r(i_r) - pipe.r_ad) * (mesh.r(i_r) - pipe.r_int)) / ((pipe.r_ad - (pipe.r_ad + pipe.r_int) / 2) * (pipe.r_int - (pipe.r_ad + pipe.r_int) / 2)));%2*salt.u_m*(1-(mesh.r(i_r)/pipe.r_int)^2);
        end

        % Initialize temperature at z = 0
  
        % Solve temperature profile along the axial direction
        z_span = mesh.z;
        [z, T] = ode15s(@(z,T) dTdz(z, T, mesh, pipe, salt, air, const), z_span, T_0, options);
        for i_z = 1:mesh.n_z
            salt.T_m(i_z, i_L_z) = sum((mesh.A_z(1,1:end-1).*salt.u_z).* T(i_z,1:end-1))/(pipe.A_c*salt.u_m);
        end

    end

figure(5);
hold on
plot(pipe.L_z_var, salt.T_m(end-1,:) - 273.15, '-o');
xlabel('Reciever Height, L [m]','FontSize',16);
ylabel('Mean Fluid Outlet Temperature [\circC]','FontSize',16);
legend({'$\dot{m}$ = 40 kg/s','$\dot{m}$ = 67 kg/s','$\dot{m}$ = 85 kg/s'}, 'Interpreter', 'latex', 'Location', 'best');
title('Parametric Study: Effect of Reciever Height on Outlet Mean Fluid Temperature');
end
%% Tube Size Parametric Study with q flux

%Standard
air.q_flux_var = [313990, 445000, 708650];
salt.mdot = 67;
    for i_qflux = 1:length(air.q_flux_var)
        
        air.q_flux = air.q_flux_var(i_qflux);
%Variable
pipe.L_z_var = 5:1:50;	% length of pipe [m]
pipe.r_int = 7.89835;	% internal radius of pipe [m]
pipe.r_ad = 7.87335; % radius of internal pipe
pipe.r_ext_var = [7.899, 7.9, 8];	% external radius of pipe [m]

for i_L_z = 1:length(pipe.L_z_var)
    pipe.L_z = pipe.L_z_var(i_L_z);
    pipe.A_c = pi*(pipe.r_int^2-pipe.r_ad^2);     % cross-sectional area for flow [m^2]
    pipe.rho = 9000;
    pipe.cp = 200;
    mesh.dr = (pipe.r_int-pipe.r_ad)/100;               % mesh radial discertization [m]
    mesh.dz = pipe.L_z/100;                % mesh axial discertization [m]
    mesh.r = (pipe.r_ad:mesh.dr:pipe.r_int);      % mesh radial locations [m]
    mesh.r = [mesh.r, pipe.r_ext];
    mesh.z = (0.0:mesh.dz:pipe.L_z)';   % mesh axial locations [m]
    mesh.n_r = length(mesh.r);               % number of radial nodes [--]
    mesh.n_z = length(mesh.z);             % number of axial nodes [--]
    %
    mesh.A_r = zeros(1,mesh.n_r);    % areas for radial conduction on the outer boundary of cell i_r [m^2]
    mesh.A_z = zeros(1,mesh.n_r);   % cross-sectional area for the z-direction flow at cell i_r [m^2]
    salt.u_z = zeros(1,mesh.n_r-1);      % axial velocity radial profile for fully developed flow.
    %
    mesh.A_z(1,1) =  pi*((mesh.r(1)+mesh.dr/2)^2 - (mesh.r(1)-mesh.dr/2)^2) ;  
    mesh.A_r(1,1) =  2*pi*mesh.dz*(mesh.r(1)+mesh.dr/2);  
    salt.u_z(1,1) =  0;%2*salt.u_m;
    for i_r = 2:mesh.n_r-2
    
      mesh.A_z(1,i_r) =  pi*((mesh.r(i_r)+mesh.dr/2)^2-(mesh.r(i_r)-mesh.dr/2)^2);  
      mesh.A_r(1,i_r) =  2*pi*mesh.dz*(mesh.r(i_r)+mesh.dr/2); 
      salt.u_z(1,i_r) =  1.5*salt.u_m * (1 - ((2*mesh.r(i_r)-(pipe.r_ad+pipe.r_int))/(pipe.r_int-pipe.r_ad)).^2); %2 * salt.u_m * (1 - ((mesh.r(i_r) - pipe.r_ad) * (mesh.r(i_r) - pipe.r_int)) / ((pipe.r_ad - (pipe.r_ad + pipe.r_int) / 2) * (pipe.r_int - (pipe.r_ad + pipe.r_int) / 2)));%2*salt.u_m*(1-(mesh.r(i_r)/pipe.r_int)^2);
    end
    mesh.A_z(1,mesh.n_r-1) =  pi*((mesh.r(mesh.n_r-1)+mesh.dr/2)^2 - (mesh.r(mesh.n_r-1)-mesh.dr/2)^2) ;   
    mesh.A_r(1,mesh.n_r-1) =  2*pi*mesh.dz*(mesh.r(mesh.n_r-1)+mesh.dr/2); 
    salt.u_z(1,mesh.n_r-1) =  0;
    
    %Adding another node for the outer surface 
    mesh.A_z(1,mesh.n_r) =  pi*((mesh.r(end)+mesh.dr/2)^2 - (mesh.r(end)-mesh.dr/2)^2) ;   
    mesh.A_r(1,mesh.n_r) =  2*pi*mesh.dz*(mesh.r(end)+mesh.dr/2); 

        salt.u_m = salt.mdot/(salt.rho*pipe.A_c); % mean velocity through pipe [m/s]
        for i_r = 2:mesh.n_r-2
                salt.u_z(1,i_r) =  1.5*salt.u_m * (1 - ((2*mesh.r(i_r)-(pipe.r_ad+pipe.r_int))/(pipe.r_int-pipe.r_ad)).^2); %2 * salt.u_m * (1 - ((mesh.r(i_r) - pipe.r_ad) * (mesh.r(i_r) - pipe.r_int)) / ((pipe.r_ad - (pipe.r_ad + pipe.r_int) / 2) * (pipe.r_int - (pipe.r_ad + pipe.r_int) / 2)));%2*salt.u_m*(1-(mesh.r(i_r)/pipe.r_int)^2);
        end

        % Initialize temperature at z = 0
  
        % Solve temperature profile along the axial direction
        z_span = mesh.z;
        [z, T] = ode15s(@(z,T) dTdz(z, T, mesh, pipe, salt, air, const), z_span, T_0, options);
        for i_z = 1:mesh.n_z
            salt.T_m(i_z, i_L_z) = sum((mesh.A_z(1,1:end-1).*salt.u_z).* T(i_z,1:end-1))/(pipe.A_c*salt.u_m);
        end

    end

figure(6);
hold on
plot(pipe.L_z_var, salt.T_m(end-1,:) - 273.15, '-o');
xlabel('Reciever Height, L [m]','FontSize',16);
ylabel('Mean Fluid Outlet Temperature [\circC]','FontSize',16);
legend({'qflux = 313990 [W/m^2]','qflux = 445000 [W/m^2]','qflux = 708650 [W/m^2]'}, 'Interpreter', 'latex', 'Location', 'best');title('Parametric Study: Effect of Reciever Height on Outlet Mean Fluid Temperature');
    end

%% Parametric Study on Wall Thickness

%Standard
air.q_flux = 445000;
salt.mdot = 67;

%Variable
pipe.L_z = 30;	% length of pipe [m]
pipe.r_int = 7.89835;	% internal radius of pipe [m]
pipe.r_ad = 7.87335; % radius of internal pipe
pipe.r_ext_var = [7.899, 7.9:0.01:9];	% external radius of pipe [m]

for i_r_ext = 1:length(pipe.r_ext_var)
    pipe.L_z = pipe.L_z;
    pipe.A_c = pi*(pipe.r_int^2-pipe.r_ad^2);     % cross-sectional area for flow [m^2]
    pipe.rho = 9000;
    pipe.cp = 200;
    mesh.dr = (pipe.r_int-pipe.r_ad)/100;               % mesh radial discertization [m]
    mesh.dz = pipe.L_z/100;                % mesh axial discertization [m]
    mesh.r = (pipe.r_ad:mesh.dr:pipe.r_int);      % mesh radial locations [m]
    mesh.r = [mesh.r, pipe.r_ext_var(i_r_ext)];
    mesh.z = (0.0:mesh.dz:pipe.L_z)';   % mesh axial locations [m]
    mesh.n_r = length(mesh.r);               % number of radial nodes [--]
    mesh.n_z = length(mesh.z);             % number of axial nodes [--]
    %
    mesh.A_r = zeros(1,mesh.n_r);    % areas for radial conduction on the outer boundary of cell i_r [m^2]
    mesh.A_z = zeros(1,mesh.n_r);   % cross-sectional area for the z-direction flow at cell i_r [m^2]
    salt.u_z = zeros(1,mesh.n_r-1);      % axial velocity radial profile for fully developed flow.
    %
    mesh.A_z(1,1) =  pi*((mesh.r(1)+mesh.dr/2)^2 - (mesh.r(1)-mesh.dr/2)^2) ;  
    mesh.A_r(1,1) =  2*pi*mesh.dz*(mesh.r(1)+mesh.dr/2);  
    salt.u_z(1,1) =  0;%2*salt.u_m;
    for i_r = 2:mesh.n_r-2
    
      mesh.A_z(1,i_r) =  pi*((mesh.r(i_r)+mesh.dr/2)^2-(mesh.r(i_r)-mesh.dr/2)^2);  
      mesh.A_r(1,i_r) =  2*pi*mesh.dz*(mesh.r(i_r)+mesh.dr/2); 
      salt.u_z(1,i_r) =  1.5*salt.u_m * (1 - ((2*mesh.r(i_r)-(pipe.r_ad+pipe.r_int))/(pipe.r_int-pipe.r_ad)).^2); %2 * salt.u_m * (1 - ((mesh.r(i_r) - pipe.r_ad) * (mesh.r(i_r) - pipe.r_int)) / ((pipe.r_ad - (pipe.r_ad + pipe.r_int) / 2) * (pipe.r_int - (pipe.r_ad + pipe.r_int) / 2)));%2*salt.u_m*(1-(mesh.r(i_r)/pipe.r_int)^2);
    end
    mesh.A_z(1,mesh.n_r-1) =  pi*((mesh.r(mesh.n_r-1)+mesh.dr/2)^2 - (mesh.r(mesh.n_r-1)-mesh.dr/2)^2) ;   
    mesh.A_r(1,mesh.n_r-1) =  2*pi*mesh.dz*(mesh.r(mesh.n_r-1)+mesh.dr/2); 
    salt.u_z(1,mesh.n_r-1) =  0;
    
    %Adding another node for the outer surface 
    mesh.A_z(1,mesh.n_r) =  pi*((mesh.r(end)+mesh.dr/2)^2 - (mesh.r(end)-mesh.dr/2)^2) ;   
    mesh.A_r(1,mesh.n_r) =  2*pi*mesh.dz*(mesh.r(end)+mesh.dr/2); 

        salt.u_m = salt.mdot/(salt.rho*pipe.A_c); % mean velocity through pipe [m/s]
        for i_r = 2:mesh.n_r-2
                salt.u_z(1,i_r) =  1.5*salt.u_m * (1 - ((2*mesh.r(i_r)-(pipe.r_ad+pipe.r_int))/(pipe.r_int-pipe.r_ad)).^2); %2 * salt.u_m * (1 - ((mesh.r(i_r) - pipe.r_ad) * (mesh.r(i_r) - pipe.r_int)) / ((pipe.r_ad - (pipe.r_ad + pipe.r_int) / 2) * (pipe.r_int - (pipe.r_ad + pipe.r_int) / 2)));%2*salt.u_m*(1-(mesh.r(i_r)/pipe.r_int)^2);
        end

        % Initialize temperature at z = 0
  
        % Solve temperature profile along the axial direction
        z_span = mesh.z;
        [z, T] = ode15s(@(z,T) dTdz(z, T, mesh, pipe, salt, air, const), z_span, T_0, options);
        for i_z = 1:mesh.n_z
            salt.T_m(i_z, i_r_ext) = sum((mesh.A_z(1,1:end-1).*salt.u_z).* T(i_z,1:end-1))/(pipe.A_c*salt.u_m);
        end

    end

figure(6);
hold on
plot(pipe.r_ext_var, salt.T_m(end-1,:) - 273.15, '-o');
xlabel('Wall Thickness, t [m]','FontSize',16);
ylabel('Mean Fluid Outlet Temperature [\circC]','FontSize',16);
%legend({'L = 10 m',}, 'Interpreter', 'latex', 'Location', 'best');
title('Parametric Study: Effect of Exterior Wall Thickness on Outlet Mean Fluid Temperature');

%% Solve for mean temperature with transience and Nusselt
%% Set up numerical mesh, radial slice cross-sectional areas, and salt axial velocity radial profile
salt.alpha = salt.k/(salt.rho*salt.c_p); % salt thermal diffusivity [m^2/s]
dt_step = 3600;     %  time interval [s]
t_max = 1*3600;   % maximum time [s]
time = 0:dt_step:t_max;   %  time [s]
mesh.time = time; %Establish time as a variable in th struc
mesh.Fo = salt.alpha*mesh.time/(pipe.r_int-pipe.r_ad); %?
% salt.T_b_mean = zeros(length(mesh.time), 1);
% salt.T_b_mean(1) = salt.T_init;  % Initial mean flow temperature

t_span = [0, t_max]; % Define the time span for integration76
% T_0 = salt.T_init * ones(mesh.n_r, 1); % Initial temperature

% Use ode15s to solve for the transient temperature profile
%[t, T] = ode15s(@(t, T) dTdt(t, T, mesh, pipe, salt, air, const), t_span, T_0, options);
% T_initial = ones(mesh.n_z, mesh.n_r) * salt.T_init;  % Assuming init everywhere initially

% Call the ODE solver (ode45)
%[t, T(mesh.n_r, mesh.n_z)] = ode45(@(t, T) dTdt(t, T, mesh, pipe, salt, air, const), t_span, T_initial);
% options_solve = optimoptions('fsolve','Algorithm','trust-region-dogleg','Display','iter',...
%     'SpecifyObjectiveGradient',false,'MaxFunEvals',200000,'MaxIter',5,'TolFun',1e-6,'StepTolerance',1e-10);
% %dTdt = @(T) dTdt(t, T, mesh, pipe, salt, air, const);
% salt.T_num = fsolve(@(T) dTdt(T, mesh, pipe, salt, air, const), T_initial, options_solve) - 273.15;
T_initial = ones(mesh.n_z, mesh.n_r) * salt.T_init;  % Assuming init everywhere initially
T_init = T_initial(:);
% options_solve = optimoptions('lsqnonlin','Display','iter',...
%     'MaxFunEvals',200000,'MaxIter',10,'TolFun',1e-6,'StepTolerance',1e-6);
% salt.T_num = lsqnonlin(@(T) dTdt(T, mesh, pipe, salt, air, const), T_initial, [], [], options_solve) - 273.15;
options = odeset('RelTol', 1e-6);
[T_transient, ~] = ode45(@(t, T) dTdt(T, mesh, pipe, salt, air, const), 0:dt_step:t_max, T_init, options);

% Reshape the solution back to the temperature grid
T_grid = reshape(T_transient(end, :), mesh.n_z, mesh.n_r);

% Plot the results
figure;
plot(mesh.time, salt.T_num);
xlabel('Time (s)');
ylabel('Temperature (K)');
legend(arrayfun(@(r) sprintf('r=%.2f', r), mesh.r, 'UniformOutput', false));
% Plot the temperature over time at a fixed z location
figure;
plot(t, T);
xlabel('Time (s)');
ylabel('Temperature at Fixed z Location');
title('Transient Temperature Profile');
% for t_i = 2:length(time)
%     t = times(t_i);
% 
%     % Update heat flux
%     q_flux_t = air.q_flux_t(t_i);
% 
%     % Update Nusselt number
%     Gz_D = salt.Re_D * salt.Pr * 2 * (pipe.r_int - pipe.r_ad) / pipe.L_z;
%     Nu_D_corr = 3.66 + 0.2362 * (Gz_D^(0.488)) * exp(-57.2 / Gz_D);
%     salt.h = Nu_D_corr * salt.k / (2 * (pipe.r_int - pipe.r_ad));
% 
%     % Compute heat input and storage term
%     q_in = air.q_flux * pipe.alpha * pi * pipe.r_ext^2;  % Solar input
%     q_out = 2 * pi * pipe.r_ext * air.h * (T_b_mean(t_idx - 1) - air.T) ...
%             + pipe.eps * const.sigma * (T_b_mean(t_idx - 1)^4 - air.T^4);  % Losses
%     q_net = q_in - q_out;
% 
%     % Update mean flow temperature using transient term
%     T_b_mean(t_idx) = T_b_mean(t_idx - 1) + ...
%                       dt * q_net / (salt.mdot * salt.c_p);
% end

%% Parametric Study on Wall Thicknes
end
%% Function that calculates Transient dTdzt for a radial profile
function [dTdt] = dTdt(T, mesh, pipe, salt, air, const)
T = reshape(T, mesh.n_z, mesh.n_r);
dTdt = zeros(mesh.n_z,mesh.n_r);
dTdt(1,:) =  0;
%air.q_flux = interp1(mesh.time, air.q_flux_t, t, 'linear', 'extrap'); % Interpolate flux based on time
m = salt.rho * mesh.dz .* mesh.A_z(1:(mesh.n_r-1));
cp = salt.c_p;
mpipe = pipe.rho * mesh.dz .* mesh.A_z(mesh.n_r);
cppipe = pipe.cp;
%z_span = mesh.z;
%T_0 = salt.T_init * ones(mesh.n_r, 1); % Initial temperature
%options = odeset('RelTol',1e-6);

% Compute dTdz using the existing steady-state calculation
%[z, T] = ode15s(@(z,T) dTdz( z, T, mesh, pipe, salt, air, const),z_span,T_0,options);

     % axial velocity radial profile for fully developed flow.
%

    % Calculate dTdt from thermal balance equation (energy conservation)
    T(1,:) = salt.T_init;
    for i_z = 2:mesh.n_z-1
        for i_r = 2:(mesh.n_r-2) % For each radial point except boundaries
            % Apply thermal balance equation: m * c_p * dT/dt = dT/dz + Q(t)
            %dTdt(i_r) = (T(i_r)) / (m(i_r) * cp);
            dTdt(i_z, i_r) = (mesh.A_z(i_r+1)*salt.rho*cp*salt.u_z(i_r+1)*(T(i_z+1,i_r)-T(i_z,i_r)) ...
        + mesh.A_z(i_r)*salt.rho*cp*salt.u_z(i_r)*(T(i_z-1,i_r)-T(i_z,i_r)) ...
        + mesh.A_r(i_r+1)*salt.k*(T(i_z,i_r+1)-T(i_z,i_r))/mesh.dr ...
        + mesh.A_r(i_r)*salt.k*(T(i_z,i_r-1)-T(i_z,i_r))/mesh.dr) ./ (m(i_r) * cp);
            %  dTdt(i_z, i_r) = (mesh.A_z(i_r+1)*salt.u_z(i_r+1)*(T(i_z+1,i_r)-T(i_z,i_r))/mesh.dz...
            % + mesh.A_z(i_r)*salt.u_z(i_r)*(T(i_z-1,i_r)-T(i_z,i_r))/mesh.dz ...
            % + mesh.A_r(i_r+1)*salt.k*(T(i_z,i_r+1)-T(i_z,i_r))/mesh.dr ...
            % + mesh.A_r(i_r)*salt.k*(T(i_z,i_r-1)-T(i_z,i_r))/mesh.dr)./(m(i_r)*cp);
            %dTdt(i_r,i_z) = mesh.A_z(i_r-1)*salt.u_z(i_r)*salt.rho*cp*
        end
        dTdt(i_z, mesh.n_r) = (pipe.alpha * air.q_flux * mesh.A_r(mesh.n_r) ...
        + pipe.k * mesh.A_z(mesh.n_r) * (T(i_z-1, mesh.n_r) - T(i_z, mesh.n_r)) ...
        - pipe.k * mesh.A_z(mesh.n_r) * (T(i_z, mesh.n_r) - T(i_z+1, mesh.n_r)) ...
        - 2 * pi * pipe.k * mesh.dz * (T(i_z, mesh.n_r) - T(i_z, mesh.n_r-1)) / log(pipe.r_ext / pipe.r_int) * mesh.A_r(mesh.n_r-1) ...
        - pipe.eps * mesh.A_r(mesh.n_r) * const.sigma * T(i_z, mesh.n_r)^4 ...
        - air.h * mesh.A_r(mesh.n_r) * (T(i_z, mesh.n_r) - air.T)) ...
        / (mpipe * cppipe);
        % dTdt(i_z,mesh.n_r) = (pipe.alpha * air.q_flux * mesh.A_r(mesh.n_r) + (pipe.k * mesh.A_z(mesh.n_r))*(T(i_z-1,mesh.n_r)-T(i_z, mesh.n_r)) - (pipe.k * mesh.A_z(mesh.n_r))*(T(i_z,mesh.n_r)-T(i_z+1, mesh.n_r)) ...
        %      - 2 * pi * pipe.k * mesh.dz * (T(i_z,mesh.n_r) - T(i_z,mesh.n_r-1)) / log(pipe.r_ext / pipe.r_int) * mesh.A_r(mesh.n_r-1) - ...
        %        pipe.eps * mesh.A_r(mesh.n_r) * const.sigma * T(:,mesh.n_r).^4 - air.h * mesh.A_r(mesh.n_r) * (T(:,mesh.n_r) - air.T)) ...
        %           / (mpipe*cppipe);
    end
    
    % Boundary conditions for the last two points (radial boundaries)
    % For node n_r-1
    dTdt(:,mesh.n_r - 1) = ((mesh.A_r(mesh.n_r - 2) * salt.k / mesh.dz) * (T(:,mesh.n_r - 2) - T(:,mesh.n_r-1)) / mesh.dr ...
                     -(mesh.A_r(mesh.n_r-1) *pipe.k / mesh.dz) * (T(:,mesh.n_r-1) - T(:,mesh.n_r)) / (pipe.r_ext-pipe.r_int)) ...
                     / (m(mesh.n_r-1)*cp);
    % For node n_r (the outer boundary)
 dTdt = dTdt(:);
end
%%  Function that calculates dTdz for an ode integrator
function [dTdz] = dTdz(z, T, mesh, pipe, salt, air, const)
dTdz = zeros(mesh.n_r,1);
dTdz(1) =  0;
for i_r = 2:(mesh.n_r-2)
    dTdz(i_r) = (mesh.A_r(i_r-1)*salt.k*((T(i_r-1)-T(i_r))/mesh.dr)...
       -(mesh.A_r(i_r)*salt.k)*((T(i_r)-T(i_r+1))/mesh.dr))...
       /(salt.rho*salt.u_z(i_r)*salt.c_p*mesh.A_z(i_r)*mesh.dz) ;
end
%
%  %  Set up derivative for boundary cell ignoring pipe axial conduction
dTdz(mesh.n_r-1) = ((mesh.A_r(mesh.n_r - 2) * salt.k / mesh.dz) * (T(mesh.n_r - 2) - T(mesh.n_r-1)) / mesh.dr ...
                     -(mesh.A_r(mesh.n_r-1) *pipe.k / mesh.dz) * (T(mesh.n_r-1) - T(mesh.n_r)) / (pipe.r_ext-pipe.r_int)) ...
                     / (mesh.A_z(mesh.n_r-1) * salt.rho * salt.u_z(mesh.n_r-2) *0.25 * salt.c_p);
dTdz(mesh.n_r) = (pipe.alpha*air.q_flux*mesh.A_r(mesh.n_r) - 2*pi*pipe.k*mesh.dz*(T(mesh.n_r)-T(mesh.n_r-1))/log(pipe.r_ext/pipe.r_int)*(mesh.A_r(mesh.n_r))-pipe.eps*mesh.A_r(mesh.n_r)*const.sigma*T(mesh.n_r)^4-air.h*mesh.A_r(mesh.n_r)*(T(mesh.n_r)-air.T))/ (mesh.A_z(mesh.n_r) *pipe.k);
end 




% function dTdz = compute_dTdz(T, mesh, pipe, salt, air, const)
%     % Calculate dTdz as you did in your steady-state function
%     dTdz = zeros(mesh.n_r, 1);
% 
%     % First cell, assuming zero gradient
%     dTdz(1) = 0;
% 
%     % Interior cells (i_r = 2 to n_r-2)
%     for i_r = 2:(mesh.n_r-2)
%         dTdz(i_r) = (mesh.A_r(i_r-1) * salt.k * ((T(i_r-1) - T(i_r)) / mesh.dr) ...
%                      - (mesh.A_r(i_r) * salt.k) * ((T(i_r) - T(i_r+1)) / mesh.dr)) ...
%                      / (salt.rho * salt.u_z(i_r) * salt.c_p * mesh.A_z(i_r) * mesh.dz);
%     end
% 
%     % Boundary condition for last two points
%     % Penultimate radial point
%     dTdz(mesh.n_r-1) = ((mesh.A_r(mesh.n_r - 2) * salt.k / mesh.dz) * (T(mesh.n_r - 2) - T(mesh.n_r - 1)) / mesh.dr ...
%                         - (mesh.A_r(mesh.n_r - 1) * pipe.k / mesh.dz) * (T(mesh.n_r - 1) - T(mesh.n_r)) / (pipe.rthe _ext - pipe.r_int)) ...
%                         / (mesh.A_z(mesh.n_r - 1) * salt.rho * salt.u_z(mesh.n_r - 2) * 0.25 * salt.c_p);
% 
%     % Outer boundary condition
%     dTdz(mesh.n_r) = (pipe.alpha * air.q_flux * mesh.A_r(mesh.n_r) - 2 * pi * pipe.k * mesh.dz * (T(mesh.n_r) - T(mesh.n_r - 1)) / log(pipe.r_ext / pipe.r_int) * mesh.A_r(mesh.n_r) ...
%                       - pipe.eps * mesh.A_r(mesh.n_r) * const.sigma * T(mesh.n_r)^4 - air.h * mesh.A_r(mesh.n_r) * (T(mesh.n_r) - air.T)) ...
%                       / (mesh.A_z(mesh.n_r) * pipe.k);
% end
