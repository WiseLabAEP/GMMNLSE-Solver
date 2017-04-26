% Tests group velocity dispersion by propagating a Gaussian pulse through
% 2*LD in SMF.
%
% Expected results: temporal broadening, peak power is reduced by a little
% more than a factor 2. See Figure 3.1 in Nonlinear Fiber Optics 3rd
% Edition, Agrawal

sim.cuda_dir_path = '../../cuda';
addpath('../../'); % MATLAB needs to know where the propagate files are located

%% Setup fiber parameters
num_modes = 1;
Aeff = 4.6263e-11;

fiber.betas = [0; 0; 24.1616/1000;]; % Dispersion coefficients, in units of ps^n/m

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

save_name = make_test_save_name('SMF_GVD', sim);

%% Setup initial conditions
N = 2^15;
tfwhm = 0.1; % ps
time_window = 20; %ps
total_energy = 0.000001;%17.04;%410.3; %nJ

L_D = (tfwhm/1.665)^2/abs(fiber.betas(3)); % dispersion length in m
fiber.L0 = 2*L_D;
sim.deltaZ = fiber.L0/100;

initial_condition = build_MMgaussian(tfwhm, time_window, total_energy, num_modes, N);

%% Run the propagation

reset(gpuDevice);
prop_output = GMMNLSE_propagate(fiber, initial_condition, sim);
save(save_name, 'prop_output', 'fiber', 'sim');
disp(prop_output.seconds);

%% Plot the results

N = size(prop_output.fields, 1);
start_I = abs(prop_output.fields(:, :, 1)).^2;
end_I = abs(prop_output.fields(:, :, end)).^2;
t = (-N/2:N/2-1)*prop_output.dt; % ps
tlim = 1;

figure();
hold on
plot(t,end_I(:, 1)/max(start_I), 'k'),axis tight, grid on
plot(t,start_I(:, 1)/max(start_I), 'k--')
ylabel('Intensity (W)')
xlabel('Time (ps)')
xlim([-tlim, tlim])
hold off