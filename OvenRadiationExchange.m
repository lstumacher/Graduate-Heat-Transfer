function[stone, oven, rad]  = MEGN571_HW6_F2024_PartD()
%
%%  MEGN571 Module 6 Fall 2024:
%       by Profs. Jackson and Porter
%


%Changes Made:

%Reduce oven height in order to increase radiation exchange between
%oven_top and the stone.

%Change virtual walls to steel walls to increase radiation exchange within
%the enclosure.

%This allows the enclosure to achieve its max temperature far faster than
%before.




close all;
clear all;
%
%% Initialize properties of the pizza stone, oven, and the external air.
sigma = 5.67e-8;                % Stefan-Boltzmann constant [W/m^2/K^4]
%
%%  MEGN571 Module 6 Assignment Fall 2024:
%       by Profs. Jackson and Porter
%
%  Solution to pizza stone radiative heat transfer in over
close all;
clear all;
%
%% Initialize properties of the pizza stone, oven, and the external air.
sigma = 5.67e-8;                % Stefan-Boltzmann constant [W/m^2/K^4]
stone.d = 0.6;                      % pizza stone diameter [m]
stone.r = 0.5*stone.d;          % pizza stone radius [m]
stone.th = 0.005;                    % pizza stone thickness [m]
%
oven.d = 0.8;                       % oven diameter [m]
oven.r = 0.5*oven.d;            % oven radius [m]
oven.L = 0.3;                       % oven height [m]
%
stone.rho = 2050;            % pizza stone density [kg/m^3]
stone.cp = 800;               % pizza stone specific heat constant [J/kg/K]
stone.k = 2.5;                  % pizza stone thermal conductivity [W/m/K]
stone.alpha = stone.k/stone.rho/stone.cp;  % pizza stone thermal diffusivity [m^2/s]
stone.eps = 0.7;             % pizza stone emissivity [--]
stone.T_init = 300;         % initial pizza stone temperature [K]
%
oven.eps = 0.8;             % oven top emissivity [--]
oven.T_top = 600;         % initial oven top temperature [K]
oven.T_wall = 300;         % outer wall oven/air temperature [K]
%
steel.eps = 0.6;              % steel ring emissivity [--]
%
% Create and index for the radiation nodes and initialize all quantities
% for radiation calculations.
stone.n_rad = 3;                % number of radiation nodes on stone
rad.ind = struct('stone', (1:stone.n_rad)', 'oven_top', stone.n_rad+1, ...
    'steel_ring', stone.n_rad+2, 'oven_wall', stone.n_rad+3, ...
    'tot', stone.n_rad+3);
%
%% Part a)  Solve for areas and view factors
[rad.F, rad.A] = ViewFactors(rad.ind, stone, oven);
%
table_names = {'stone_1','stone_2','stone_3','oven_top','steel_ring','oven_wall'};
Table_F = array2table(rad.F,'VariableNames',table_names,'RowNames',table_names);
Table_F
%
%% Part b) Set up radiosity equations and solve for values at initial temperatures
%
% Set up stone ring volumes and masses
stone.A = rad.A(1:stone.n_rad);
stone.Vol = stone.A*stone.th;
stone.m = stone.Vol*stone.rho;
%
% Allocate memory of radiation vectors for calculating radiosities and heat fluxes
rad.eps = zeros(rad.ind.tot,1);
rad.eps = [stone.eps; stone.eps; stone.eps; oven.eps; steel.eps; steel.eps];
rad.q = zeros(rad.ind.tot,1);
rad.T = [stone.T_init; stone.T_init; stone.T_init; oven.T_top; stone.T_init; oven.T_wall];
E_b = zeros(rad.ind.tot,1);
J_guess = zeros(rad.ind.tot,1);
%
% Initialize surfaces emissivities, initial temperatures, and blackbodies in rad structure

% Initalizing the blackbody Emmisive Power via E = sigma*T^4
E_b = sigma*[stone.T_init;stone.T_init;stone.T_init;oven.T_top; oven.T_wall; oven.T_wall].^4; % ???
%Guess values for sigma that are similar to the emissive power
J_guess = sigma*[stone.T_init;stone.T_init;stone.T_init; oven.T_wall; oven.T_wall; oven.T_wall].^4; % ???
%
% Solve for initial radiositites at time = 0 and initial temperatures
options = optimoptions('fsolve','Algorithm','trust-region-dogleg','Display','none');
rad.J_0 = fsolve(@(J) J_calc(J, E_b, rad), J_guess, options);
%
rad.J_0
%
%% Part c) Integrate transient temperature equations to get stone T(t)
%
%  Set up time span for 60 min
time_f = 60*60;                             % time of simulation [s]
t_span = linspace(0, time_f, 60);  % time span for reporting [s]
T_0 = rad.T(rad.ind.stone);           % initial temperature for transient integration [K]
%
%  Perform time integration based on radiative heat exchange (ignoring conduction)
[t, T] = ode45(@(t,T) stone_dTdt(t, T, sigma, rad, stone), t_span, T_0);
%
%  Make plot of temperature with respect to time for all of the pizze stone
r = linspace(0, stone.r, stone.n_rad);
plot(t/60, T() - 273.15, 'LineWidth',2) 
xlabel('Time in oven [min]')
ylabel('Stone temperature [deg. C]')
legend('N=1','Box','Off','Location','northwest')
legend(append('r = ' ,num2str(r(1)*100),' cm'), ...
    append('r = ' ,num2str(r(2)*100),' cm'), ...
    append('r = ' ,num2str(r(3)*100),' cm'), ...
    'Box','Off','Location','southeast')
set(gca,'FontSize',18)
end
%
%% Function to update stone temperature in stoneT for current time step
function [dTdt] = stone_dTdt(t, T, sigma, rad, stone)
%
% Update black body radiation terms
E_b = sigma*[T; rad.T(rad.ind.oven_top); rad.T(rad.ind.steel_ring); rad.T(rad.ind.oven_wall)].^4;
%
% Solve for radiosities
options = optimoptions('fsolve','Algorithm','trust-region-dogleg','Display','none');
J = fsolve(@(J) J_calc(J, E_b, rad), E_b, options);
%
%  Calculate heat fluxes for surfaces given radiosities
[q] =  Jtoq(J, E_b, rad);
%
%  Calculate temperature derivatives for pizza stone rigns

dTdt = (-q(rad.ind.stone))./(stone.m*stone.cp); %stone.eps*sigma*stone.A.*(600^4-T.^4)./(stone.m*stone.cp);%(q(rad.ind.stone))./(stone.m*stone.cp);  % ???;.
%dTdt = dTdt(:);  % Convert to a column vector

end
%
%% Function that calculates radiosities from known heat transfer rates
function[Res] = J_calc(J, E_b, rad)
%
% Calculate sums in radiosity equations
sum_rad = zeros(rad.ind.tot,1);

%
% ???
%
% Calculate a residual to solve for radiosity for stone rings

Res(rad.ind.stone(1)) = (E_b(1) - J(rad.ind.stone(1)))/((1-rad.eps(1))/(rad.eps(1)*rad.A(rad.ind.stone(1))))...  
    - rad.A(rad.ind.stone(1))*rad.F(rad.ind.stone(1),rad.ind.oven_top)*(J(rad.ind.stone(1))-J(rad.ind.oven_top)) ...
    - rad.A(rad.ind.stone(1))*rad.F(rad.ind.stone(1),rad.ind.oven_wall)*(J(rad.ind.stone(1))-J(rad.ind.oven_wall)); % ???

Res(rad.ind.stone(2)) = (E_b(2) - J(rad.ind.stone(2)))/((1-rad.eps(2))/(rad.eps(2)*rad.A(rad.ind.stone(2))))...  
    - rad.A(rad.ind.stone(2))*rad.F(rad.ind.stone(2),rad.ind.oven_top)*(J(rad.ind.stone(2))-J(rad.ind.oven_top)) ...
    - rad.A(rad.ind.stone(2))*rad.F(rad.ind.stone(2),rad.ind.oven_wall)*(J(rad.ind.stone(2))-J(rad.ind.oven_wall)); % ???

Res(rad.ind.stone(3)) = (E_b(3) - J(rad.ind.stone(3)))/((1-rad.eps(3))/(rad.eps(3)*rad.A(rad.ind.stone(3))))... 
    - rad.A(rad.ind.stone(3))*rad.F(rad.ind.stone(3),rad.ind.oven_top)*(J(rad.ind.stone(3))-J(rad.ind.oven_top)) ...
    - rad.A(rad.ind.stone(3))*rad.F(rad.ind.stone(3),rad.ind.oven_wall)*(J(rad.ind.stone(3))-J(rad.ind.oven_wall)); % ???

%
% Calculate residuals to solve for radiosity for metal ring, oven_top, and outer oven wall

Res(rad.ind.steel_ring) =  0 ...
    - rad.A(rad.ind.steel_ring)*rad.F(rad.ind.steel_ring,rad.ind.oven_top)*(J(rad.ind.steel_ring)-J(rad.ind.oven_top))...
    - rad.A(rad.ind.steel_ring)*rad.F(rad.ind.steel_ring,rad.ind.oven_wall)*(J(rad.ind.steel_ring)-J(rad.ind.oven_wall)); % ???

Res(rad.ind.oven_top) = (E_b(4) - J(rad.ind.oven_top))/((1-rad.eps(4))/(rad.eps(4)*rad.A(rad.ind.oven_top)))...
    - rad.A(rad.ind.oven_top)*rad.F(rad.ind.oven_top,rad.ind.stone(1))*(J(rad.ind.oven_top)-J(rad.ind.stone(1))) ...
    - rad.A(rad.ind.oven_top)*rad.F(rad.ind.oven_top,rad.ind.stone(2))*(J(rad.ind.oven_top)-J(rad.ind.stone(2))) ...
    - rad.A(rad.ind.oven_top)*rad.F(rad.ind.oven_top,rad.ind.stone(3))*(J(rad.ind.oven_top)-J(rad.ind.stone(3))) ... 
    - rad.A(rad.ind.oven_top)*rad.F(rad.ind.oven_top,rad.ind.steel_ring)*(J(rad.ind.oven_top)-J(rad.ind.steel_ring)) ...
    - rad.A(rad.ind.oven_top)*rad.F(rad.ind.oven_top,rad.ind.oven_wall)*(J(rad.ind.oven_top)-J(rad.ind.oven_wall)); % ???

Res(rad.ind.oven_wall) = 0 ... %(E_b(6) - J(rad.ind.oven_wall))/((1-rad.eps(6))/(rad.eps(6)*rad.A(rad.ind.oven_top))) ...
    - rad.A(rad.ind.oven_wall)*rad.F(rad.ind.oven_wall,rad.ind.stone(1))*(J(rad.ind.oven_wall)-J(rad.ind.stone(1))) ...
    - rad.A(rad.ind.oven_wall)*rad.F(rad.ind.oven_wall,rad.ind.stone(2))*(J(rad.ind.oven_wall)-J(rad.ind.stone(2))) ...
    - rad.A(rad.ind.oven_wall)*rad.F(rad.ind.oven_wall,rad.ind.stone(3))*(J(rad.ind.oven_wall)-J(rad.ind.stone(3))) ...
    - rad.A(rad.ind.oven_wall)*rad.F(rad.ind.oven_wall,rad.ind.steel_ring)*(J(rad.ind.oven_wall)-J(rad.ind.steel_ring)) ...
    - rad.A(rad.ind.oven_wall)*rad.F(rad.ind.oven_wall,rad.ind.oven_top)*(J(rad.ind.oven_wall)-J(rad.ind.oven_top)); % ???
end
%
%% Function to calculate surface temperatures or heat fluxes from radiosities
function[q] = Jtoq(J, E_b, rad)
%
% Initialize heat fluxes
q = zeros(rad.ind.tot,1);

% Calculate heat fluxes for pizza stone central circle and outer rings
q(rad.ind.stone(1)) = (E_b(rad.ind.stone(1)) - J(rad.ind.stone(1))) / (1-rad.eps(1)) * (rad.eps(1)*rad.A(rad.ind.stone(1)));
q(rad.ind.stone(2)) = (E_b(rad.ind.stone(2)) - J(rad.ind.stone(2))) / (1-rad.eps(2)) * (rad.eps(2)*rad.A(rad.ind.stone(2)));
q(rad.ind.stone(3)) = (E_b(rad.ind.stone(3)) - J(rad.ind.stone(3))) / (1-rad.eps(3)) * (rad.eps(3)*rad.A(rad.ind.stone(3)));

% Calculate heat flux from oven top
q(rad.ind.oven_top) = (E_b(rad.ind.oven_top) - J(rad.ind.oven_top)) / (1-rad.eps(4)) * (rad.eps(4)*rad.A(rad.ind.oven_top));

% Calculate heat flux for the steel ring
q(rad.ind.steel_ring) = (E_b(rad.ind.steel_ring) - J(rad.ind.steel_ring)) / (1-rad.eps(5)) * (rad.eps(5)*rad.A(rad.ind.steel_ring));

% Calculate summation for heat flux on the oven wall
q(rad.ind.oven_wall) = rad.A(rad.ind.oven_wall) * ...
    (rad.F(rad.ind.oven_wall, rad.ind.stone(1)) * (J(rad.ind.oven_wall) - J(rad.ind.stone(1))) + ...
     rad.F(rad.ind.oven_wall, rad.ind.stone(2)) * (J(rad.ind.oven_wall) - J(rad.ind.stone(2))) + ...
     rad.F(rad.ind.oven_wall, rad.ind.stone(3)) * (J(rad.ind.oven_wall) - J(rad.ind.stone(3))) + ...
     rad.F(rad.ind.oven_wall, rad.ind.oven_top) * (J(rad.ind.oven_wall) - J(rad.ind.oven_top)) + ...
     rad.F(rad.ind.oven_wall, rad.ind.steel_ring) * (J(rad.ind.oven_wall) - J(rad.ind.steel_ring)));

% sigma = 5.67e-8;
% % Initialize heat fluxes
% q = zeros(rad.ind.tot,1);
% %
% % Calculate heat fluxes for pizza stone central circle and outer rings
% 
% q(rad.ind.stone(1)) = (E_b(rad.ind.stone(1))  - J(rad.ind.stone(1)))/((1-rad.eps(1))/(rad.eps(1)*rad.A(rad.ind.stone(1))));
% q(rad.ind.stone(2)) = (E_b(rad.ind.stone(2)) - J(rad.ind.stone(2)))/((1-rad.eps(2))/(rad.eps(2)*rad.A(rad.ind.stone(2))));
% q(rad.ind.stone(2)) = (E_b(rad.ind.stone(3)) - J(rad.ind.stone(3)))/((1-rad.eps(3))/(rad.eps(3)*rad.A(rad.ind.stone(3))));
% 
% 
% %
% % Calcalate heat flux from oven top, temperature for adiabatic ring, and
% % heat flux for external wall
% 
% q(rad.ind.oven_top) = (E_b(rad.ind.oven_top) - J(rad.ind.oven_top))/((1-rad.eps(4)))/(rad.eps(4)*rad.A(rad.ind.oven_top)); % 
% 
% 
% q(rad.ind.steel_ring) = (E_b(rad.ind.steel_ring) - J(rad.ind.steel_ring)) / (1-rad.eps(5)) * (rad.eps(5)*rad.A(rad.ind.steel_ring));; % ???
% %T(rad.ind.steel_ring) = (J(rad.ind.steel_ring)/sigma)^0.25;
% 
% % Calculate suumation for heat flux calculation to get oven wall calcultion
% %
% %  ???
% %
% q(rad.ind.oven_wall) = rad.A(rad.ind.oven_wall)*rad.F(rad.ind.oven_wall,rad.ind.stone(1))*(J(rad.ind.oven_wall)-J(rad.ind.stone(1))) ...
%     + rad.A(rad.ind.oven_wall)*rad.F(rad.ind.oven_wall,rad.ind.stone(2))*(J(rad.ind.oven_wall)-J(rad.ind.stone(2))) ...
%     + rad.A(rad.ind.oven_wall)*rad.F(rad.ind.oven_wall,rad.ind.stone(3))*(J(rad.ind.oven_wall)-J(rad.ind.stone(3))) ...
%     + rad.A(rad.ind.oven_wall)*rad.F(rad.ind.oven_wall,rad.ind.steel_ring)*(J(rad.ind.oven_wall)-J(rad.ind.steel_ring)) ...
%     + rad.A(rad.ind.oven_wall)*rad.F(rad.ind.oven_wall,rad.ind.oven_top)*(J(rad.ind.oven_wall)-J(rad.ind.oven_top));
end
%
%% Function to calculate view factors pizza oven
function[F, A] = ViewFactors(ind, stone, oven)
%
%  Allocate memory for areas and view factors for radiation nodes
A = zeros(ind.tot,1);
F = zeros(ind.tot,ind.tot);
%
% Set inner and outer radii for node rings (see C-47 configuration factor)
dr = stone.r/(stone.n_rad);
r_1 = linspace(0, stone.r-dr, stone.n_rad);
r_2 = linspace(dr, stone.r, stone.n_rad);
%
% Gind areas for different radiation nodes
A(ind.stone) = pi*(r_2.^2 - r_1.^2); % area of rings on pizza stone [m^2]
A(ind.oven_top) = pi*oven.r^2;    % area of oven ceiling [m^2]
A(ind.steel_ring) = pi*(oven.r^2 - stone.r^2); % area of steel ring [m^2]
A(ind.oven_wall) = 2*pi*oven.r*oven.L;    % area of oven wall [m^2]
%
% Find view factors from oven top to pizza stone rings
% ???
R1 = dr/oven.L;
Rtop = oven.r/oven.L;
X41 = 1+(1+R1^2)/Rtop^2;
F(ind.oven_top,ind.stone(1)) = 0.5*(X41-((X41^2)-4*(R1/Rtop)^2)^0.5); % ???
%
% Find view factors from oven top to pizza stone rings
% ???
for i_rad = 2:ind.stone(end)
    Htop = oven.L/oven.r;
    R2 = dr*(i_rad-1)/oven.r;
    R3 = dr*(i_rad)/oven.r;
    F(ind.oven_top,ind.stone(i_rad)) = 0.5*(R3^2 - R2^2 - ((1+ R3^2 + Htop^2)^2 - 4* R3^2)^0.5 + ((1+R2^2+Htop^2)^2 - 4* R2^2)^0.5); % ???
end
%
% Fnd view factor from oven top to steel ring
% ???
R5 = oven.r/oven.L;
Rtop = oven.r/oven.L;
X45 = 1+(1+R5^2)/Rtop^2;
F(ind.oven_top,ind.steel_ring) = 0.5*(X45-((X45^2)-4*(R5/Rtop)^2)^0.5) - F(ind.oven_top,ind.stone(1))-F(ind.oven_top,ind.stone(2))-F(ind.oven_top,ind.stone(3)); % ???
%
% Use reciprocity for view factors from oven ceiling to stone and steel surfaces (rad.ind.oven_top)
% and use summation for view factors from 1:rad.ind.steel_ring to wall
for i_rad = 1:ind.tot-1
    F(i_rad,ind.oven_top) = A(ind.oven_top)*F(ind.oven_top,i_rad)/A(i_rad);  % ???
    if i_rad == 4
        F(ind.oven_top,ind.oven_wall) = 1 -  F(ind.oven_top,ind.stone(1)) - F(ind.oven_top,ind.stone(2)) - F(ind.oven_top,ind.stone(3))- F(ind.oven_top,ind.oven_top) - F(ind.oven_top,ind.steel_ring);
    else 
    F(i_rad,ind.oven_wall) = 1 - F(i_rad, ind.oven_top);
    end
    F(ind.oven_wall, i_rad) = A(i_rad)*F(i_rad,ind.oven_wall)/A(ind.oven_wall);  % ???
end
F(ind.oven_wall,ind.oven_wall) = 1 - F(ind.oven_wall,1)- F(ind.oven_wall,2) - F(ind.oven_wall,3) - F(ind.oven_wall,4) - F(ind.oven_wall,5); % ???
end