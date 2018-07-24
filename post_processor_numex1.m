clc
close all
clear all

pre_define_the_path_for_folders

% This routine loads the mat file generated at the end of simulation
% and prints/plots the output data

disp('-----------------------------------------------------')
disp('Post-processing...')

% ------------------------------------------------------------------------
%                     Printing the R+jL matrices
% -------------------------------------------------------------------------

load('results_numex1_straight_conductor/data_R_jL_mat.mat');
disp('-----------------------------------------------------')
disp('R+jL matrices ::: ')

for freq_no=1:num_freq
    disp(['Frequency = ',num2str(freq_all(freq_no))])
    for kk=1:num_ports
        disp([num2str(R_jL_mat(kk,:,freq_no))])
    end
end
cd results_numex1_straight_conductor
fasthenry_results_straight_conductor
cd ..

FigHandle = figure;
% delay visibility as long as possible, as this speeds up visualization
set(FigHandle, 'Visible', 'off');

set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
%set(FigHandle, 'Position', [100, 100, 1280, 1024]);
set(FigHandle, 'Position', [100, 100, 1280, 900]);
subplot(2,1,1)
set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
h=loglog(freq_all,1e3*real(squeeze(R_jL_mat)),'b-'); set(h,'LineWidth',2);
hold on
h=loglog(fH_data(:,1),1e3*real(fH_data(:,2)),'ro'); set(h,'LineWidth',2);
legend('VoxHenry','FastHenry');
axis tight;grid on;xlabel('Frequency (Hz)');ylabel('Resistance (m{\Omega})');
ylim(1e3*[0.005 0.021]); set(gca,'YTick',1e3*[0.005 0.01 0.015 0.02]);
set(gca,'XTick',[10^0 10^2 10^4 10^6 10^8 10^10])
set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
set(gca,'LineWidth',1); grid minor;

%FigHandle = figure;
set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
%set(FigHandle, 'Position', [100, 100, 1280, 1024]);
%set(FigHandle, 'Position', [100, 100, 1280, 900]);
subplot(2,1,2)
set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
h=loglog(freq_all,1e12*imag(squeeze(R_jL_mat)),'b-'); set(h,'LineWidth',2);
hold on
h=loglog(fH_data(:,1),1e12*imag(fH_data(:,2)),'ro'); set(h,'LineWidth',2);
legend('VoxHenry','FastHenry');
axis tight;grid on;xlabel('Frequency (Hz)');ylabel('Inductance (pH)');
if(exist ("OCTAVE_VERSION", "builtin") > 0)
    % need to work-around Octave limit/bug with low ylim values
    ylim(gca, 'auto');
else
    ylim(1e12*[0.95e-11 1.06e-11]); set(gca,'YTick',1e12*[0.95e-11 0.975e-11 1.0e-11 1.025e-11 1.05e-11]);
end
set(gca,'XTick',[10^0 10^2 10^4 10^6 10^8 10^10])
set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');
set(gca,'LineWidth',1); grid minor;

set(FigHandle, 'Visible', 'on');
%refresh;
drawnow;

print('results_numex1_straight_conductor/res_ind_voxhenry',  '-dpng', '-r300')
print('Results/res_ind_voxhenry',  '-dpng', '-r300')


disp('L2 norm error for resistance :::')
sqrt((sum(abs(real(squeeze(R_jL_mat))-real(fH_data(:,2))).^2)/sum(abs(real(fH_data(:,2))).^2)))

disp('L2 norm error for inductance :::')
sqrt(sum(abs(imag(squeeze(R_jL_mat))-imag(fH_data(:,2))).^2)/sum(abs(imag(fH_data(:,2))).^2))


disp('-----------------------------------------------------')

% ------------------------------------------------------------------------
%                     Printing CPU Times
% -------------------------------------------------------------------------

load('results_numex1_straight_conductor/data_CPU_timings.mat');
disp('-----------------------------------------------------')
disp('CPU times ::: ')
disp(['Time for generating Ae matrix ::: ',num2str(sim_CPU_pre(1))])
disp(['Time for generating circulant tensors + RHS vector ::: ',num2str(sim_CPU_pre(2))])
tot_prep=sum(sim_CPU_pre);
disp(['Total Time for preparing LSE data ::: ',num2str(tot_prep)])

tot_sol=zeros(num_freq,1);
for freq_no=1:num_freq
    disp(['For frequency = ',num2str(freq_all(freq_no))])
    for port_no=1:num_ports
        if (port_no == 1)
            disp(['Time for FFT of circulant ::: ',num2str(sim_CPU_lse(freq_no,port_no,1))])
            tot_sol(freq_no)=tot_sol(freq_no)+sim_CPU_lse(freq_no,port_no,1);
            disp(['Time for generating sparse precon ::: ',num2str(sim_CPU_lse(freq_no,port_no,2))])
            tot_sol(freq_no)=tot_sol(freq_no)+sim_CPU_lse(freq_no,port_no,2);
        end
        disp(['Time for iterative solution for port #',num2str(port_no),' ::: ',num2str(sim_CPU_lse(freq_no,port_no,3))])
        tot_sol(freq_no)=tot_sol(freq_no)+sim_CPU_lse(freq_no,port_no,3);
    end
    disp(['Total Time for solving freq pnt ',num2str(freq_no),'::: ',num2str(tot_sol(freq_no))])
end

disp('Summary ::: ')
disp(['Total Time for preparing LSE data ::: ',num2str(tot_prep)])
disp(['Total Time for solving LSE ::: ',num2str(sum(tot_sol))])
disp(['Total Time for simulation ::: ',num2str(tot_prep+sum(tot_sol))])
disp('-----------------------------------------------------')

% ------------------------------------------------------------------------
%                     Plotting Current Distribution
% -------------------------------------------------------------------------
% select plotting option - check the subroutines below for more options
% option 1-> total currents on 3D structure, - no cut selection required
% option 2-> current coefficients on the voxels - select plane and cut
% option 3-> currents on the nodes via imagesc - select plane and cut
% option 4-> currents on the nodes via quiver - select plane and cut
% option 5-> currents on the structure w/directions via quiver3 - no cut selection required
% voxels (on a selected cut), 3-> currents on nodes w/scalar values (on a selected cut)
load('results_numex1_straight_conductor/data_curr_plot.mat')
plot_option=1;
[L,M,N] = size(Mc);

disp('-----------------------------------------------------')
disp(['Plotting Current Distribution...'])

% if any of plot option 2,3,4 is selected, define plane and cut
slct_plane='xy'; %'xz'; 'yz';
if (plot_option == 2 || plot_option == 3)
    % 1) use the following for plot option 2 and 3
    slct_cut=round(N/2);% round(M/2); round(L/2);
elseif (plot_option == 4)
    % 2) use the following for plot option 4 - we need coordinate of the cut
    slct_cut=squeeze(r(1,1,N,3)); % z-coordinate of cut % squeeze(r(round(L/2),1,1,1)); % x-coordinate of cut; squeeze(r(1,round(M/2),1,2)); % y-coordinate of cut
end

if (plot_option == 2)
    % sort current coefficients on voxels
    [Jx_currs_grid,Jy_currs_grid,Jz_currs_grid,J2d_currs_grid,J3d_currs_grid,cmin,cmax]=post_obtain_curr_coefs_on_grid(x,Mc);
elseif (plot_option > 2)
    % obtain currents on nodes
    [nodes_w_currs_x_aligned,nodes_w_currs_y_aligned,nodes_w_currs_z_aligned]=post_obtain_currs_on_nodes(x,Ae_only_leaving,Ae_only_entering_bndry,r,Mc,dx);
end

switch plot_option
    case 1
        % 1) Plot currents on structure
        % Plot total currents as one scalar on each voxel
        plot_currs_on_3D_structure(x,Ae_only_leaving,r,Mc,dx)
    case 2
        % 2) Plot current coefficients obtained via iterative solution
        plot_curr_coefs_on_grid(slct_plane,slct_cut,r,Jx_currs_grid,Jy_currs_grid,Jz_currs_grid,J2d_currs_grid,J3d_currs_grid,cmin,cmax);
    case 3
        % 3) Plot currents with scalar values via imagesc
        plot_curr_on_nodes(slct_plane,slct_cut,dx,nodes_w_currs_x_aligned,nodes_w_currs_y_aligned,nodes_w_currs_z_aligned)
    case 4
        % 4) Plot currents on cuts w/ directions via quiver
        plot_curr_on_nodes_quiver(slct_plane,slct_cut,nodes_w_currs_x_aligned,nodes_w_currs_y_aligned,nodes_w_currs_z_aligned)
    case 5
        % 5) Plot currents on the structure w/directions via quiver3
        plot_curr_on_nodes_quiver3(nodes_w_currs_x_aligned,nodes_w_currs_y_aligned,nodes_w_currs_z_aligned)
    otherwise
        disp('No current plotting!')
end

print('results_numex1_straight_conductor/curr_dist',  '-dpng', '-r300')
print('Results/curr_dist',  '-dpng', '-r300')

disp(['Done... Plotting Current Distribution'])
disp('-----------------------------------------------------')

disp('Done... Post-processing')
disp('-----------------------------------------------------')

