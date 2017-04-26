function output = build_MMgaussian(tfwhm, time_window, total_energy, num_modes, N)
% build_MMgaussian  Build a multimode gaussian from parameters
% tfwhm - full witdth at half maximum, in ps
% time_window - full time window width in ps
% total_energy - total energy of the pulse, in all modes, in nJ
% num_modes - number of modes
% N - number of time grid points

t0 = tfwhm/1.665; % ps
dt = time_window/N; % ps
t = (-N/2:N/2-1)*dt; % ps

% electric field envelope in W^0.5
time_profile = sqrt(total_energy/num_modes/(t0*sqrt(pi))*10^3)*exp(-t.^2/(2*t0^2));

field = zeros(N, num_modes);
for idx = 1:num_modes
    field(:, idx) = time_profile;
end

% output a struct
output.fields = field;
output.dt = dt;

end