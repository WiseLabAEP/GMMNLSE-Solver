Fibers

This folder is a good place to store fibers that are built with the combination of solve_for_modes, calc_dispersion, and calc_SRSK_tensors. You can find an example of how to build a fiber in Examples/build_fiber_example, which builds a MM fiber that supports 6 modes at 1030 nm.

In general this set of steps solves for the modes at a number of wavelengths so the dispersion coefficients can be generated, however once the coefficients have been calculated one no longer needs the modes at all wavelengths. To save space, therefore, the fibers that have been included here by default only keep the modes at 1030 nm. To see the full output with the modes at all wavelengths, run Examples/build_fiber_example.