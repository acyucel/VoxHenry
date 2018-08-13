%function print_out_inputs_generate_consts
if (num_freq > 1)
    if (issorted(freq) == 0) % not sorted
        freq=sort(freq)
    end
    freq_all = freq;
    freq = freq(1); % currently do everything for the lowest freq
elseif (num_freq == 1)
    freq_all = freq;
end
% generate EM constants
EMconstants;
disp('------------------------------------------------------------------')
disp('VoxHenry: Inductance Extraction Simulator for Voxelized Geometries')
disp('                                                                  ')
disp('          by Abdulkadir C. Yucel and Jacob K. White (MIT)         ')
disp('                                                                  ')
disp('                 Empowered by modules/ideas from                  ')
disp('   Athanasios G. Polimeridis and Ioannis P. Georgakis (Skoltech)  ')
disp('     Hakan Bagci (KAUST), Enrico Di Lorenzo (FastFieldSolvers)    ')
disp('------------------------------------------------------------------')
disp('Inputs for simulation :::')
if(num_freq == 1)
    disp(['Frequency = ',num2str(freq)]);
    disp(['Conductivity of conductors = ',num2str(se)]);
    tmp_er=er - 1j*se/(eo*omega);
    disp(['Relative epsr of conductors (Re/Im) = ',num2str(real(tmp_er)),' / ',num2str(-imag(tmp_er))])
    skin_depth=1/sqrt(pi*freq*4*pi*1e-7*se);
    disp(['Discretization should be at order of (skin depth) ', num2str(skin_depth)])
else
    disp(['Frequencies = ']);
    for kk=1:num_freq
        fprintf('%10.5e ', freq_all(kk))
    end
    disp('  ')
    disp(['Conductivity of conductors = ',num2str(se)]);
    tmp_er=er - 1j*se/(eo*(2*pi*freq_all(1)));
    disp(['Relative epsr of conductors (Re/Im) - at min freq = ',num2str(real(tmp_er)),' / ',num2str(-imag(tmp_er))])
    tmp_er=er - 1j*se/(eo*(2*pi*freq_all(end)));
    disp(['Relative epsr of conductors (Re/Im) - at max freq = ',num2str(real(tmp_er)),' / ',num2str(-imag(tmp_er))])
    skin_depth=1/sqrt(pi*freq_all(1)*4*pi*1e-7*se);
    disp(['Skin depth at min freq = ', num2str(skin_depth)])
    skin_depth=1/sqrt(pi*freq_all(end)*4*pi*1e-7*se);
    disp(['Skin depth at max freq = ', num2str(skin_depth)])
end

disp(['Voxel Size = ', num2str(Res)])
disp(['Number of Ports = ', num2str(num_ports)])
disp(['Tolerance for iterative solver = ', num2str(tol)])
disp(['# of maximum inner / outer GMRES iterations = ', num2str(inner_it),' / ',num2str(outer_it)])

if (fl_check_ports == 1 || fl_check_geo == 1 || fl_check_ports == 1)
    disp('***************** No Simulation will be performed *****************')
    disp('Just for checking the geometry, computational domain, or port nodes')
    disp('Turn of fl_check_ports, fl_check_geo, fl_check_ports flags for simulation')
end

if (fl_check_geo == 1)
    disp('Offline mode: The geometry will be plotted ')
end
if (fl_check_domain == 1)
    disp('Offline mode: The computational domain will be plotted ')
end
if (fl_check_ports == 1)
    disp('Offline mode: The port nodes will be plotted ')
end

% store CPU times for the simulations
sim_CPU_pre=zeros(1,2);
% last entry: 1-> Ae time, 2-> circulant tensor + rhs generation time
sim_CPU_lse=zeros(num_freq,num_ports,3); 
% last entry: 1-> FFT and dataprep, 2-> sparse precon, 3-> iter sol

% check whether the freq_curr_plot is in freq list 
if(num_freq == 1) % automatically assign single frequency value
    freq_curr_plot = freq;
else % otherwise search in frequency list
    fl_found = 0;
    for kk=1:num_freq
        if (abs(freq_all(kk) - freq_curr_plot) < 1e-12)
            fl_found = 1;
        end
        if (fl_found == 1)
            break
        end
    end
    
    if (fl_found == 0) % if it is not in the list, assign the value to the last one
        freq_curr_plot = freq(end);
    end
end

if (plot_option >= 1 && plot_option <= 1 )
    disp(['Currents will be plotted at the frequency ', num2str(freq_curr_plot)])
end
