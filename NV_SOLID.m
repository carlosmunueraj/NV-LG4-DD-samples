%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Final sensor simulations                          %
%                       Carlos Munuera Javaloy                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% -things to load: bz, Bz, Omeff, OmLG, DLG, oscpoints, reps
load('NewNuclearSignals/Solid_DD_RK4_150kHz_OU.mat')

% Constants
%--------------------------------------------------------------------------
ge = 28.024e9;    % Hz/T
gh = 42.577478e6; % Nuclear gyromagnetic ratio (Hz/Tesla)
 % Nuclear gyromagnetic ratio (Hz/Tesla)
hbar = 1.054571628e-34;  % Hbar
T = 300;
kB = 1.38064852e-23;   % Julios * Kelvin-1 
D = 2.87e9;
therm = 2*pi*hbar*gh*Bz/(kB*T);
alpha = 55*pi/180;
density = 5.2e28; %Detectable proton density in ethanol
geom = 4.1; %Geometric integral
mu0 = 4e-7*pi;  % N/A^2

dt = 1/(Omeff*oscpoints);

Om = 40e6;
tpi = 1/(2*Om);

% Crosstalk
OmNV = OmLG * ge / gh;
wLnuc = -gh * Bz;
pp = D-ge*Bz;
wd_plus = wLnuc-DLG; % Driving frequency for the positive axes
wd_minus = wLnuc+DLG; % Driving frequency for the negative axes
d_plus = (pp-wd_plus); % Detuning for positive axes
d_minus = (pp-wd_minus); % Detuning for negative axes

%REPS5 REMOVED
% Parameters
%--------------------------------------------------------------------------
padding = 0*1e5;

% Simulation stuff
rp = pi/2-acos(sqrt(3)*cos(alpha)/(sqrt(2+cos(2*alpha))));
pa = round(rp / (2*pi) * oscpoints);
pb = pa+round(oscpoints/2);
pT = oscpoints;

% Operators
%--------------------------------------------------------------------------
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];
sz = [1 0; 0 -1];

sphi = @(phi) sx*cos(phi)+sy*sin(phi);

H0 = @(bz, ct) -bz*ge*sz/2+ct*sz/2;

rho0 = [1 1; 1 1]/2;

Ux = expm(-1i*pi*sx/2);

tht_plane = acos(1/sqrt(2*cos(alpha)^2+1));
Uf = expm(1i*pi/2*sx/2);

% Inner things
%--------------------------------------------------------------------------
seq = [pa, pb, 2*pT-pb, 2*pT-pa, 2*pT+pa, 2*pT+pb];
seq_phi = [0, 0, 0, 0, 0, 0];
obs = {sz, sz, -sz, -sz};
Cross = {OmNV^2/(4*d_plus), OmNV^2/(4*d_minus), OmNV^2/(4*d_minus), OmNV^2/(4*d_plus)};

t_vec = (1:reps)*4/Omeff-1/Omeff;

% Loop
%--------------------------------------------------------------------------
measurements = zeros(1, reps);

ind_seq = 0;
for r = 1:reps
    ind_phi = 1;
    ind0 = (r-1)*4*oscpoints+1;

    phi = 0;

    rho = rho0;
    for l = 1:3*oscpoints
        ct = Cross{ceil(l/oscpoints)};
        %ct=0;
        U = expm(-1i*2*pi*H0(bz(ind0 + l - 1), ct)*dt);
        
        if ismember(l, seq)
            U = expm(-1i*2*pi*(H0(bz(ind0 + l - 1), ct)+Om*sphi(phi)/2)*tpi)*expm(-1i*2*pi*H0(bz(ind0 + l - 1), ct)*(dt-tpi));
        end
        rho = U*rho*U';
    end
    rho = Uf*rho*Uf';
    measurements(r) = real(trace(rho*sz));
end

measurements = measurements-mean(measurements);

%% Fourier
%--------------------------------------------------------------------------
[fspectrum, frecs] = fourier([measurements, zeros(1, padding)], 4/Omeff);

%% Plots
%--------------------------------------------------------------------------
figure(1)
coef = 5*12*ge*hbar^2*gh^2*mu0*density*Bz*geom/(32*Omeff*pi^2*kB*T)*sin(atan(sqrt(2)));
coef = (2*pi)^3*coef;
CS = [3.66, 1.19]*1e-6;
delta_star = CS*Bz*gh*sqrt(1+2*cos(alpha)^2)/3;
[meas] = analytical(coef, delta_star, t_vec);

subplot(2, 1, 1);
hold on
plot(t_vec, measurements, 'LineWidth', 2)
plot(t_vec, (meas), 'LineWidth', 2)
xlim([0,0.5])

subplot(2, 1, 2);
hold on
plot(frecs, imag(fspectrum), 'LineWidth', 2)

% Dressing and checks
CS = [3.66, 1.19]*1e-6;
xline(CS*Bz*gh*sqrt(1+2*cos(alpha)^2)/3, 'LineWidth', 2)
xlim([0, 500])
%ylim([0, 7e-5])

% Auxiliar functions
%--------------------------------------------------------------------------
function [fspectrum, frecs] = fourier(X, dt)
    L = length(X);
    Y = fftshift(fft(X));
    
    frecs = ((-L/2):((L/2)-1))/(L*dt);
    fspectrum = Y/L;
end

% Values used in AERIS: err_amp=0.0024, tau=0.001
function [eps] = OU(eps0, tau, err_amp, dt)
    eps = eps0*exp(-dt/tau)+err_amp*randn;
end

function [meas] = analytical(coef, delta_star, t_vec)

    meas = .6*coef*cos(2*pi*delta_star(1)*t_vec-pi/2)...
         + .4*coef*cos(2*pi*delta_star(2)*t_vec-pi/2);
end