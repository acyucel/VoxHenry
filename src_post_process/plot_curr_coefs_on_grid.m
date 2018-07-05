function plot_curr_coefs_on_grid(slct_plane,slct_cut,r,Jx_currs_grid,Jy_currs_grid,Jz_currs_grid,J2d_currs_grid,J3d_currs_grid,cmin,cmax)


fl_ind_or_phys_len=1; % if 1 axes are physical length, otherwise index along that dimension
fl_clbar_to_maxmin=0; % if 1 colorbar limits are set to max and min currents in [Jx,Jy,Jz] 
fl_diag_currs=1;

if (fl_ind_or_phys_len == 1)
    x_var_grid=squeeze(r(:,1,1,1));
    y_var_grid=squeeze(r(1,:,1,2));
    z_var_grid=squeeze(r(1,1,:,3));
else
    x_var_grid=1:1:length(squeeze(r(:,1,1,1)));
    y_var_grid=1:1:length(squeeze(r(1,:,1,2)));
    z_var_grid=1:1:length(squeeze(r(1,1,:,3)));
end

fl_slct_part=3; %1-real,2-imag,3-abs
if (fl_slct_part == 1) % real part of data for plotting
    Jx_currs_grid=real(Jx_currs_grid);Jy_currs_grid=real(Jy_currs_grid);
    Jz_currs_grid=real(Jz_currs_grid);J2d_currs_grid=real(J2d_currs_grid);
    J3d_currs_grid=real(J3d_currs_grid);
    str_Jx='{\it{Re}}(J_x) ';str_Jy='{\it{Re}}(J_y) ';str_Jz='{\it{Re}}(J_z) ';
    str_J2d='{\it{Re}}(J_{2D}) '; str_J3d='{\it{Re}}(J_{3D}) ';
elseif (fl_slct_part == 2) % imag part of data for plotting
    Jx_currs_grid=imag(Jx_currs_grid);Jy_currs_grid=imag(Jy_currs_grid);
    Jz_currs_grid=imag(Jz_currs_grid);J2d_currs_grid=imag(J2d_currs_grid);
    J3d_currs_grid=imag(J3d_currs_grid);
    str_Jx='{\it{Im}}(J_x) ';str_Jy='{\it{Im}}(J_y) ';str_Jz='{\it{Im}}(J_z) ';
    str_J2d='{\it{Im}}(J_{2D}) '; str_J3d='{\it{Im}}(J_{3D}) ';
elseif (fl_slct_part == 3)% abs value of data for plotting
    Jx_currs_grid=abs(Jx_currs_grid);Jy_currs_grid=abs(Jy_currs_grid);
    Jz_currs_grid=abs(Jz_currs_grid);J2d_currs_grid=abs(J2d_currs_grid);
    J3d_currs_grid=abs(J3d_currs_grid);
    str_Jx='|J_x| ';str_Jy='|J_y| ';str_Jz='|J_z| ';
    str_J2d='|J_{2D}| '; str_J3d='|J_{3D}| ';
end




if (slct_plane == 'xz')
    hor_axis='z';
    vert_axis='x';
    
    figure; set(gca,'FontName','Times New Roman');
    subplot(2,2,1); set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    imagesc(z_var_grid,x_var_grid,squeeze(Jx_currs_grid(:,slct_cut,:))); title(str_Jx);colorbar; ylabel(vert_axis); xlabel(hor_axis); set(gca,'FontSize',24); if(fl_clbar_to_maxmin == 1);caxis([cmin cmax]); end;
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    subplot(2,2,2); set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    imagesc(z_var_grid,x_var_grid,squeeze(Jy_currs_grid(:,slct_cut,:))); title(str_Jy);colorbar; ylabel(vert_axis); xlabel(hor_axis); set(gca,'FontSize',24);if(fl_clbar_to_maxmin == 1);caxis([cmin cmax]); end;
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    subplot(2,2,3); set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    imagesc(z_var_grid,x_var_grid,squeeze(Jz_currs_grid(:,slct_cut,:))); title(str_Jz);colorbar; ylabel(vert_axis); xlabel(hor_axis); set(gca,'FontSize',24);if(fl_clbar_to_maxmin == 1);caxis([cmin cmax]); end;
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    if (fl_diag_currs == 1)
        figure;
        subplot(2,2,1); set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
        imagesc(z_var_grid,x_var_grid,squeeze(J2d_currs_grid(:,slct_cut,:))); title(str_J2d);colorbar; ylabel(vert_axis); xlabel(hor_axis); set(gca,'FontSize',24);
        set(gca,'FontSize',24); set(gca,'FontName','Times New Roman'); if(fl_clbar_to_maxmin == 1);caxis([cmin cmax]); end;
        subplot(2,2,2); set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
        imagesc(z_var_grid,x_var_grid,squeeze(J3d_currs_grid(:,slct_cut,:))); title(str_J3d);colorbar; ylabel(vert_axis); xlabel(hor_axis); set(gca,'FontSize',24);
        set(gca,'FontSize',24); set(gca,'FontName','Times New Roman'); if(fl_clbar_to_maxmin == 1);caxis([cmin cmax]); end;
    end
elseif (slct_plane == 'xy')  % cut parallel to xy
    
    hor_axis='y';
    vert_axis='x';
    
    figure;
    subplot(2,2,1); set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    imagesc(y_var_grid,x_var_grid,squeeze(Jx_currs_grid(:,:,slct_cut))); title(str_Jx);colorbar; ylabel(vert_axis); xlabel(hor_axis); set(gca,'FontSize',24); if(fl_clbar_to_maxmin == 1);caxis([cmin cmax]); end;
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    subplot(2,2,2); set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    imagesc(y_var_grid,x_var_grid,squeeze(Jy_currs_grid(:,:,slct_cut))); title(str_Jy);colorbar; ylabel(vert_axis); xlabel(hor_axis); set(gca,'FontSize',24);if(fl_clbar_to_maxmin == 1);caxis([cmin cmax]); end;
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    subplot(2,2,3); set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    imagesc(y_var_grid,x_var_grid,squeeze(Jz_currs_grid(:,:,slct_cut))); title(str_Jz);colorbar; ylabel(vert_axis); xlabel(hor_axis); set(gca,'FontSize',24);if(fl_clbar_to_maxmin == 1);caxis([cmin cmax]); end;
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    
    if (fl_diag_currs == 1)
        figure;
        subplot(2,2,1); set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
        imagesc(y_var_grid,x_var_grid,squeeze(J2d_currs_grid(:,:,slct_cut))); title(str_J2d);colorbar; ylabel(vert_axis); xlabel(hor_axis); set(gca,'FontSize',24);
        set(gca,'FontSize',24); set(gca,'FontName','Times New Roman'); if(fl_clbar_to_maxmin == 1);caxis([cmin cmax]); end;
        subplot(2,2,2); set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
        imagesc(y_var_grid,x_var_grid,squeeze(J3d_currs_grid(:,:,slct_cut))); title(str_J3d);colorbar; ylabel(vert_axis); xlabel(hor_axis); set(gca,'FontSize',24);
        set(gca,'FontSize',24); set(gca,'FontName','Times New Roman'); if(fl_clbar_to_maxmin == 1);caxis([cmin cmax]); end;
    end
elseif (slct_plane == 'yz') % cut parallel to yz
    
    hor_axis='z';
    vert_axis='y';
    
    figure;
    subplot(2,2,1); set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    imagesc(z_var_grid,y_var_grid,squeeze(Jx_currs_grid(slct_cut,:,:))); title(str_Jx);colorbar; ylabel(vert_axis); xlabel(hor_axis); set(gca,'FontSize',24); if(fl_clbar_to_maxmin == 1);caxis([cmin cmax]); end;
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    subplot(2,2,2); set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    imagesc(z_var_grid,y_var_grid,squeeze(Jy_currs_grid(slct_cut,:,:))); title(str_Jy);colorbar; ylabel(vert_axis); xlabel(hor_axis); set(gca,'FontSize',24); if(fl_clbar_to_maxmin == 1);caxis([cmin cmax]); end;
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    subplot(2,2,3); set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    imagesc(z_var_grid,y_var_grid,squeeze(Jz_currs_grid(slct_cut,:,:))); title(str_Jz);colorbar; ylabel(vert_axis); xlabel(hor_axis); set(gca,'FontSize',24); if(fl_clbar_to_maxmin == 1);caxis([cmin cmax]); end;
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    if (fl_diag_currs == 1)
        figure;
        subplot(2,2,1); set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
        imagesc(z_var_grid,y_var_grid,squeeze(J2d_currs_grid(slct_cut,:,:))); title(str_J2d);colorbar; ylabel(vert_axis); xlabel(hor_axis); set(gca,'FontSize',24);
        set(gca,'FontSize',24); set(gca,'FontName','Times New Roman'); if(fl_clbar_to_maxmin == 1);caxis([cmin cmax]); end;
        subplot(2,2,2); set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
        imagesc(z_var_grid,y_var_grid,squeeze(J3d_currs_grid(slct_cut,:,:))); title(str_J3d);colorbar; ylabel(vert_axis); xlabel(hor_axis); set(gca,'FontSize',24);
        set(gca,'FontSize',24); set(gca,'FontName','Times New Roman'); if(fl_clbar_to_maxmin == 1);caxis([cmin cmax]); end;
    end
end