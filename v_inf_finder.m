%% V-inf Finder
% by Parham Khodadi
clear; close all; clc;

%% Conversions
yr_to_sec = 31556952;       % Years to seconds
AU_to_m = 149597870700;     % Astronomical units to meters

%% Givens
muSun   = 1.32712440042e20;     % [m^3/s^2]

%% Inputs
rp_AU   = input('What is the perihelion distance? [AU]: ');         % Perihelion radius [AU]
rf_AU   = input('What is the target distance? [AU]: ');             % target radius [AU]
T_yr    = input('How long should the missions take? [years]: ');    % flight time [years]

%% Calculations
rp      = rp_AU * AU_to_m;      % Perihelion radius [m]
T       = T_yr * yr_to_sec;     % flight time [s]
rf      = rf_AU * AU_to_m;      % target radius [m]
a_from_vinf   = @(vinf) muSun./(vinf.^2);
e_from_vinf   = @(vinf) 1 + rp./(a_from_vinf(vinf));  % e > 1
F_from_vinf   = @(vinf) acosh( (rf./(a_from_vinf(vinf)) + 1)./e_from_vinf(vinf) );
tof_from_vinf = @(vinf) sqrt((a_from_vinf(vinf)).^3/muSun) .* ...
                        ( e_from_vinf(vinf).*sinh(F_from_vinf(vinf)) - F_from_vinf(vinf) );

% Find v_inf such that TOF = T aka TOF - T = 0
f = @(vinf) tof_from_vinf(vinf) - T;
vinf_guess = rf/T;                      % v_avg
vinf = fzero(f, vinf_guess);            % heliocentric hyperbolic excess [m/s]

fprintf('v_inf_sun  = %.4f km/s\n', vinf/1e3);

v_p = sqrt(vinf^2 + 2*muSun/rp);
fprintf('v_perihelion  = %.4f km/s\n', v_p/1e3);