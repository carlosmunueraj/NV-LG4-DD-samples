%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Sample simulation with RK4 and OU                     %
%                       Carlos Munuera Javaloy                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

folder_name = "NewNuclearSignals/";
file_name = "SOLID_DD_RK4_100kHz_OU.mat";

% Constants
%--------------------------------------------------------------------------
gh = 42.577478e6; % Nuclear gyromagnetic ratio (Hz/Tesla)
Bz = 2.1;
hbar = 1.054571628e-34;  % Hbar
mu0 = 4e-7*pi;  % N/A^2
T = 300;
kB = 1.38064852e-23;   % Julios * Kelvin-1 
therm = 2*pi*hbar*gh*Bz/(kB*T);
alpha = 55*pi/180;
density = 5.2e28; %Detectable proton density in ethanol
geom = 4.1; %Geometric integral

% Parameters
%--------------------------------------------------------------------------
% System
nuc_coords = ([-0.178/sqrt(3),      0,        0;
                0.178/(2*sqrt(3)),  0.178/2,  0;
                0.178/(2*sqrt(3)), -0.178/2,  0;
                0.178/sqrt(3),      0,        sqrt(0.248^2-0.178^2/3);
               -0.178/(2*sqrt(3)),  0.178/2,  sqrt(0.248^2-0.178^2/3)]*1e-9)';

N = size(nuc_coords, 2); % Number of nuclei
euler = [0, pi/5, pi/3]; % Rotate the molecule about ZYZ
CS = [3.66, 3.66, 3.66, 1.19, 1.19]; % Scalar chemical shifts

% Drivings
OmLG = 100e3;
DLG = OmLG / sqrt(2);
Omeff = sqrt(OmLG^2 + DLG^2);

% Relaxation
T2 = 0.2;

% Simulation
t_tot = 0.5; % Total simulation time
oscpoints = 100; % Points per oscillation

dtLG = 1 / (Omeff * oscpoints);
reps = round(t_tot*Omeff/4);
pointsLG = 4*reps*oscpoints;

% Noise
Navg = 20;

tau_OU = 1e-3;
percent = 0.0024;
c = 2*(OmLG*percent)^2/tau_OU;

% Operators
%--------------------------------------------------------------------------
[Ix, Iy, Iz, Id, Itotx, Itoty, Itotz] = operators(N);

% Prepare rho
rho_i = [exp(2*pi*hbar*gh*Bz/(2*kB*T)) 0; 0 exp(-2*pi*hbar*gh*Bz/(2*kB*T))]...
       /(exp(-2*pi*hbar*gh*Bz/(2*kB*T))+exp(2*pi*hbar*gh*Bz/(2*kB*T)));
rho0 = rho_i;
for k = 1:(N-1)
    rho0 = kron(rho0, rho_i);
end

[Hcs, Hdd] = Hmol(nuc_coords, CS, euler, Ix, Iy, Iz, Bz, gh, mu0, hbar);
H0 = Hcs+Hdd;

Hx = @(noise) (OmLG+noise) * Itotx;
Hy = @(noise) (OmLG+noise) * Itoty;
Hd = DLG * Itotz;

%Accomodate the error study
%Hamiltonians
HA =@(noise) H0 + Hx(noise) * sin(alpha) + Hy(noise) * cos(alpha) + Hd;
HmA =@(noise) H0 - Hx(noise) * sin(alpha) - Hy(noise) * cos(alpha) - Hd;
HmB =@(noise) H0 + Hx(noise) * sin(alpha) - Hy(noise) * cos(alpha) - Hd;
HB =@(noise) H0 - Hx(noise) * sin(alpha) + Hy(noise) * cos(alpha) + Hd;

% Master equation
ME = @(rho, H) -1i*2*pi*(H*rho-rho*H)+...
     1/(2*T2)*(4*(Iz(:, :, 1)*rho*Iz(:, :, 1)+Iz(:, :, 2)*rho*Iz(:, :, 2)+...
     Iz(:, :, 3)*rho*Iz(:, :, 3)+Iz(:, :, 4)*rho*Iz(:, :, 4)+...
     Iz(:, :, 5)*rho*Iz(:, :, 5))-N*rho);

% Evolution
%--------------------------------------------------------------------------
tpi0 = (pi/2-acos(1/sqrt(2*cos(alpha)^2+1)))/(2*pi*OmLG);

pi0_points = round(tpi0/dtLG);
dt_pi0 = tpi0/pi0_points;

rhoy = rho0;
for k = 1:pi0_points
    rhoy = RK4(rhoy, @(rho) ME(rho, Hx(0)+H0), dt_pi0);
end

tic
spectrumLG = zeros(1, pointsLG);
parfor m=1:Navg
    rho = rhoy;
    temparray = zeros(1, pointsLG);
    
    noise = OU(0, 1, tau_OU, c);

    ind = 1;
    for r = 1:reps
        noise = OU(noise, 4/OmLG, tau_OU, c);
        MEs = {@(rho) ME(rho, HA(noise)), @(rho) ME(rho, HmA(noise)),...
               @(rho) ME(rho, HmB(noise)), @(rho) ME(rho, HB(noise))};
    
        for l = 1:4
            MEfun = MEs{l};
    
            for k = 1:oscpoints
                rho = RK4(rho, MEfun, dtLG);
                temparray(ind) = real(trace(rho*Itotz/N));
    
                ind = ind + 1;
            end
        end
    end
    spectrumLG = spectrumLG+temparray;
end
toc
%%
%bz = spectrumLG/Navg*1.3e-9*kB*T*2/(N*pi*hbar*gh);%Careful we changed the
%observable from Itotz to Itotz/N
bz = 2*pi*gh*hbar*mu0*density*geom/(4*pi)*spectrumLG/Navg;
save(strcat(folder_name, file_name), "bz", "Bz", "Omeff", "OmLG", "DLG", "oscpoints", "reps");

%% Plot
%--------------------------------------------------------------------------
hold on
plot(bz)

%% Auxiliary functions
%--------------------------------------------------------------------------
function [y_f] = RK4(y_i, f, h)
    k1 = f(y_i);
    k2 = f(y_i+h*k1/2);
    k3 = f(y_i+h*k2/2);
    k4 = f(y_i+h*k3);

    y_f = y_i+h/6*(k1+2*k2+2*k3+k4);
end

function [noise] = OU(noise, dt, tau_OU, c)
    noise = noise*exp(-dt/tau_OU)+sqrt(c*tau_OU/2*(1-exp(-2*dt/tau_OU)))*randn;
end