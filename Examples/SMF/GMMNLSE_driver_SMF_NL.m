% Tests nonlinear propagation without dispersion.
%
% Expected results: Spectral broadening. See Figure 4.2 in Nonlinear Fiber Optics 3rd
% Edition, Agrawal

sim.cuda_dir_path = '../../cuda';
addpath('../../'); % MATLAB needs to know where the propagate files are located

%% Setup fiber parameters
num_modes = 1;
Aeff = 4.6263e-11;

fiber.betas = [0; 0; 0; 0]; % Dispersion coefficients, in units of ps^n/m

SR = ones(1, 1, 1, 1);
SR(1, 1, 1, 1) = 1/Aeff;
fiber.SR = SR;


%% Setup simulation parameters
c = 2.99792458e-4; %speed of ligth m/ps
lambda = 1030e-9; % m

sim.f0=c/lambda; % central pulse frequency (THz)
sim.fr = 0;%0.18;
sim.sw = 0;
sim.M = 10;
sim.n_tot_max = 20;
sim.n_tot_min = 2;
sim.tol = 5*10^-4;
sim.save_period = 0; % Just set it to be the fiber length
sim.SK_factor = 1;
sim.check_nan = 1;
sim.verbose = 1;
if ~isfield(sim, 'defaults_set') || sim.defaults_set == 0
    sim.single_yes = 1;
    sim.gpu_yes = 1;
    sim.mpa_yes = 1;
end

save_name = make_test_save_name('SMF_NL', sim);

%% Setup initial conditions
N = 2^15;
tfwhm = 0.1; % ps
time_window = 20; %ps
total_energy = 1;%17.04;%410.3; %nJ

w0 = 2*pi*sim.f0; % angular frequency (THz)
n2 = 2.3*10^-20; % m^2 W^-1
gamma = n2*w0/(Aeff*c); % W^-1 m

initial_condition = build_MMgaussian(tfwhm, time_window, total_energy, num_modes, N);

L_NL = 1/(gamma*max(abs(initial_condition.fields).^2)); % dispersion length in m
fiber.L0 = 5*pi*L_NL;
sim.deltaZ = fiber.L0/1000;
sim.save_period = fiber.L0/10;

%% Run the propagation

reset(gpuDevice);
prop_output = GMMNLSE_propagate(fiber, initial_condition, sim);
save(save_name, 'prop_output', 'fiber', 'sim');
disp(prop_output.seconds);

%% Plot the results

N = size(prop_output.fields, 1);
I_freq0 = abs(ifftshift(ifft(prop_output.fields(:, :, 1)))).^2; % z = 0
I_freqpi = abs(ifftshift(ifft(prop_output.fields(:, :, 3)))).^2; % z = pi*L_NL
I_freq3p5pi = abs(ifftshift(ifft(prop_output.fields(:, :, 8)))).^2; % z = 3.5*pi*L_NL
f = sim.f0+(-N/2:N/2-1)/(prop_output.dt*N); % ps
flim = 50;

figure();
hold on
plot(f,I_freq0, 'k'),axis tight, grid on
plot(f,I_freqpi, 'b'),axis tight, grid on
plot(f,I_freq3p5pi, 'r'),axis tight, grid on
ylabel('Intensity (a.u.)')
xlabel('Frequency (THz)')
xlim([sim.f0-flim, sim.f0+flim])
hold off