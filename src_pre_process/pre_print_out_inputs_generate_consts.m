disp('Inputs for simulation :::')
% take as reference 'se' the maximum conductivity of any conductor
se = max(max(max(sigma_e)));
% 'er' is instead zero
er = 0;

disp(['Domain dimension: ', num2str(L), ' x ', num2str(M), ' x ', num2str(N), ' = ', num2str(L*M*N), ' voxels']);
disp(['Number of non-empty voxels: ', num2str(size(idxS,1))]);
disp(['Voxel Size = ', num2str(dx)])
disp(['Number of Ports = ', num2str(num_ports)])

if(num_freq == 1)
    disp(['Frequency = ',num2str(freq)]);
    disp(['Conductivity of conductors = ',num2str(se)]);
    tmp_er=er - 1j*se/(eo*omega);
    disp(['Relative epsr of conductors (Re/Im) = ',num2str(real(tmp_er)),' / ',num2str(-imag(tmp_er))])
    if se > 0.0
      skin_depth=1/sqrt(pi*freq*4*pi*1e-7*se);
      disp(['Skin depth at given frequency = ', num2str(skin_depth)])
    end
else
    fprintf('Frequencies = ');
    for kk=1:num_freq
        fprintf('%10.5e ', freq_all(kk))
    end
    disp('  ')
    disp(['Conductivity of conductors = ',num2str(se)]);
    tmp_er=er - 1j*se/(eo*(2*pi*freq_all(1)));
    disp(['Relative epsr of conductors (Re/Im) - at min freq = ',num2str(real(tmp_er)),' / ',num2str(-imag(tmp_er))])
    tmp_er=er - 1j*se/(eo*(2*pi*freq_all(end)));
    disp(['Relative epsr of conductors (Re/Im) - at max freq = ',num2str(real(tmp_er)),' / ',num2str(-imag(tmp_er))])
    if se > 0.0
      skin_depth=1/sqrt(pi*freq_all(1)*4*pi*1e-7*se);
      disp(['Skin depth at min freq = ', num2str(skin_depth)])
      skin_depth=1/sqrt(pi*freq_all(end)*4*pi*1e-7*se);
      disp(['Skin depth at max freq = ', num2str(skin_depth)])
    end
end
if ~isempty(lambdaL)
  lL = max(max(max(lambdaL)));
  disp(['London penetration depth (max) = ', num2str(lL)])
end    

disp(['Tolerance for iterative solver = ', num2str(tol)])
disp(['# of maximum inner / outer GMRES iterations = ', num2str(inner_it),' / ',num2str(outer_it)])

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

if (plot_option == 1)
    disp(['Currents will be plotted at the frequency ', num2str(freq_curr_plot)])
end
