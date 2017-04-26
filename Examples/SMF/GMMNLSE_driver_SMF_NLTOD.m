% Tests combined effects of nonlinearity and TOD
%
% Expected results: Oscillating tail in time, broadened spectrum with more
% energy towards the red. See Figure 4.14 in Nonlinear Fiber Optics 3rd
% Edition, Agrawal

sim.cuda_dir_path = '../../cuda';
addpath('../../'); % MATLAB needs to know where the propagate files are located

%% Setup fiber parameters
num_modes = 1;
Aeff = 4.6263e-11;

fiber.betas = [0; 0; 0; 1000/1000^2]; % Dispersion coefficients, in units of ps^n/m

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

save_name = make_test_save_name('SMF_NLTOD', sim);

%% Setup initial conditions
N = 2^15;
tfwhm = 0.1; % ps
time_window = 20; %ps

N_bar = 1;
L_TOD = (tfwhm/1.665)^3/abs(fiber.betas(4)); % dispersion length in m
L_NL = L_TOD/sqrt(N_bar);

w0 = 2*pi*sim.f0; % angular frequency (THz)
n2 = 2.3*10^-20; % m^2 W^-1
gamma = n2*w0/(Aeff*c); % W^-1 m

total_energy = 1/(gamma*L_NL)*(tfwhm/1.665)*sqrt(pi)/10^3;

initial_condition = build_MMgaussian(tfwhm, time_window, total_energy, num_modes, N);

fiber.L0 = 5*L_TOD;
sim.deltaZ = fiber.L0/1000;

%% Run the propagation

reset(gpuDevice);
prop_output = GMMNLSE_propagate(fiber, initial_condition, sim);
save(save_name, 'prop_output', 'fiber', 'sim');
disp(prop_output.seconds);

%% Plot the results

N = size(prop_output.fields, 1);
I_time = abs(prop_output.fields(:, :, end).^2);
I_freq = abs(ifftshift(ifft(prop_output.fields(:, :, end)))).^2;
t = (-N/2:N/2-1)*(prop_output.dt);
f = sim.f0+(-N/2:N/2-1)/(prop_output.dt*N); % ps
tlim = 1;
flim = 20;

figure();
subplot(1, 2, 1);
plot(t, I_time, 'k'),axis tight, grid on
ylabel('Intensity (W)')
xlabel('Time (ps)')
xlim([-tlim, tlim])

subplot(1, 2, 2);
plot(f, I_freq, 'k'),axis tight, grid on
ylabel('Intensity (a.u.)')
xlabel('Frequency (THz)')
xlim([sim.f0-flim, sim.f0+flim])