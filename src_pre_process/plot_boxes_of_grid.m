function plot_boxes_of_grid(r,res)

fl_plot_w_voxel_cens = 0;

figure

xd = r(:,:,:,1);
yd = r(:,:,:,2);
zd = r(:,:,:,3);

if (fl_plot_w_voxel_cens == 1)
    h=plot3(xd(:), yd(:), zd(:), 'r*');
    set(h,'MarkerSize',10);
end

[L,M,N,~] = size(r); % domain size
nD = L*M*N;

num_elem_allowed = 10000;

if (nD > num_elem_allowed)% visualization of boxes is slow when nD > num_elem_allowed
    fl_method = 0; % turn off the box visualization
else
    fl_method = 2; % turn on box visualization: method (1) traditional, (2) facet based
end

fl_method = 2;

if (fl_method > 0)
    
    x_unit=[1 0 0];
    y_unit=[0 1 0];
    z_unit=[0 0 1];
    
    %if (L > 1); res=r(2,1,1,1)-r(1,1,1,1);
    %elseif (M > 1); res=r(1,2,1,2)-r(1,1,1,2)
    %elseif (N > 1); res=r(1,1,2,3)-r(1,1,1,3);
    %end
    
    half_res=0.5*res;
    
    if (fl_method == 1)
        
        box_corners=cell(size(r,1),size(r,2),size(r,3));
        for kk=1:size(r,1)
            for ll=1:size(r,2)
                for mm=1:size(r,3)
                    temp_coor=squeeze(r(kk,ll,mm,1:3))';
                    box_corners{kk,ll,mm}(1,1:3)=temp_coor-y_unit*half_res-x_unit*half_res-z_unit*half_res;
                    box_corners{kk,ll,mm}(2,1:3)=temp_coor-y_unit*half_res+x_unit*half_res-z_unit*half_res;
                    box_corners{kk,ll,mm}(3,1:3)=temp_coor+y_unit*half_res+x_unit*half_res-z_unit*half_res;
                    box_corners{kk,ll,mm}(4,1:3)=temp_coor+y_unit*half_res-x_unit*half_res-z_unit*half_res;
                    
                    box_corners{kk,ll,mm}(5,1:3)=temp_coor-y_unit*half_res-x_unit*half_res+z_unit*half_res;
                    box_corners{kk,ll,mm}(6,1:3)=temp_coor-y_unit*half_res+x_unit*half_res+z_unit*half_res;
                    box_corners{kk,ll,mm}(7,1:3)=temp_coor+y_unit*half_res+x_unit*half_res+z_unit*half_res;
                    box_corners{kk,ll,mm}(8,1:3)=temp_coor+y_unit*half_res-x_unit*half_res+z_unit*half_res;
                end
            end
        end
        
        % plotting lines between corners
        bnd_slcted_hori=['-'];
        bnd_slcted_vert=[':'];
        ln_wdth=2;
        dum_shft=4;
        clr_hori='b';
        clr_vert='b';
        
        for kk=1:size(r,1)
            for ll=1:size(r,2)
                for mm=1:size(r,3)
                    if (sum(abs(r(kk,ll,mm,1:3)))~=0)
                        line([box_corners{kk,ll,mm}(1,1),box_corners{kk,ll,mm}(2,1)],[box_corners{kk,ll,mm}(1,2),box_corners{kk,ll,mm}(2,2)],[box_corners{kk,ll,mm}(1,3),box_corners{kk,ll,mm}(2,3)] ,'LineStyle',bnd_slcted_hori,'LineWidth',ln_wdth,'Color',clr_hori);
                        line([box_corners{kk,ll,mm}(2,1),box_corners{kk,ll,mm}(3,1)],[box_corners{kk,ll,mm}(2,2),box_corners{kk,ll,mm}(3,2)],[box_corners{kk,ll,mm}(2,3),box_corners{kk,ll,mm}(3,3)] ,'LineStyle',bnd_slcted_hori,'LineWidth',ln_wdth,'Color',clr_hori);
                        line([box_corners{kk,ll,mm}(3,1),box_corners{kk,ll,mm}(4,1)],[box_corners{kk,ll,mm}(3,2),box_corners{kk,ll,mm}(4,2)],[box_corners{kk,ll,mm}(3,3),box_corners{kk,ll,mm}(4,3)] ,'LineStyle',bnd_slcted_hori,'LineWidth',ln_wdth,'Color',clr_hori);
                        line([box_corners{kk,ll,mm}(4,1),box_corners{kk,ll,mm}(1,1)],[box_corners{kk,ll,mm}(4,2),box_corners{kk,ll,mm}(1,2)],[box_corners{kk,ll,mm}(4,3),box_corners{kk,ll,mm}(1,3)] ,'LineStyle',bnd_slcted_hori,'LineWidth',ln_wdth,'Color',clr_hori);
                        
                        line([box_corners{kk,ll,mm}(dum_shft+1,1),box_corners{kk,ll,mm}(dum_shft+2,1)],[box_corners{kk,ll,mm}(dum_shft+1,2),box_corners{kk,ll,mm}(dum_shft+2,2)],[box_corners{kk,ll,mm}(dum_shft+1,3),box_corners{kk,ll,mm}(dum_shft+2,3)] ,'LineStyle',bnd_slcted_hori,'LineWidth',ln_wdth,'Color',clr_hori);
                        line([box_corners{kk,ll,mm}(dum_shft+2,1),box_corners{kk,ll,mm}(dum_shft+3,1)],[box_corners{kk,ll,mm}(dum_shft+2,2),box_corners{kk,ll,mm}(dum_shft+3,2)],[box_corners{kk,ll,mm}(dum_shft+2,3),box_corners{kk,ll,mm}(dum_shft+3,3)] ,'LineStyle',bnd_slcted_hori,'LineWidth',ln_wdth,'Color',clr_hori);
                        line([box_corners{kk,ll,mm}(dum_shft+3,1),box_corners{kk,ll,mm}(dum_shft+4,1)],[box_corners{kk,ll,mm}(dum_shft+3,2),box_corners{kk,ll,mm}(dum_shft+4,2)],[box_corners{kk,ll,mm}(dum_shft+3,3),box_corners{kk,ll,mm}(dum_shft+4,3)] ,'LineStyle',bnd_slcted_hori,'LineWidth',ln_wdth,'Color',clr_hori);
                        line([box_corners{kk,ll,mm}(dum_shft+4,1),box_corners{kk,ll,mm}(dum_shft+1,1)],[box_corners{kk,ll,mm}(dum_shft+4,2),box_corners{kk,ll,mm}(dum_shft+1,2)],[box_corners{kk,ll,mm}(dum_shft+4,3),box_corners{kk,ll,mm}(dum_shft+1,3)] ,'LineStyle',bnd_slcted_hori,'LineWidth',ln_wdth,'Color',clr_hori);
                        
                        line([box_corners{kk,ll,mm}(1,1),box_corners{kk,ll,mm}(5,1)],[box_corners{kk,ll,mm}(1,2),box_corners{kk,ll,mm}(5,2)],[box_corners{kk,ll,mm}(1,3),box_corners{kk,ll,mm}(5,3)] ,'LineStyle',bnd_slcted_vert,'LineWidth',ln_wdth,'Color',clr_vert);
                        line([box_corners{kk,ll,mm}(2,1),box_corners{kk,ll,mm}(6,1)],[box_corners{kk,ll,mm}(2,2),box_corners{kk,ll,mm}(6,2)],[box_corners{kk,ll,mm}(2,3),box_corners{kk,ll,mm}(6,3)] ,'LineStyle',bnd_slcted_vert,'LineWidth',ln_wdth,'Color',clr_vert);
                        line([box_corners{kk,ll,mm}(3,1),box_corners{kk,ll,mm}(7,1)],[box_corners{kk,ll,mm}(3,2),box_corners{kk,ll,mm}(7,2)],[box_corners{kk,ll,mm}(3,3),box_corners{kk,ll,mm}(7,3)] ,'LineStyle',bnd_slcted_vert,'LineWidth',ln_wdth,'Color',clr_vert);
                        line([box_corners{kk,ll,mm}(4,1),box_corners{kk,ll,mm}(8,1)],[box_corners{kk,ll,mm}(4,2),box_corners{kk,ll,mm}(8,2)],[box_corners{kk,ll,mm}(4,3),box_corners{kk,ll,mm}(8,3)] ,'LineStyle',bnd_slcted_vert,'LineWidth',ln_wdth,'Color',clr_vert);
                    end
                end
            end
        end
        
        axis tight;grid on;xlabel('x');ylabel('y');zlabel('z');
        
    elseif (fl_method == 2)

        for kk=1:size(r,1)
            for ll=1:size(r,2)
                for mm=1:size(r,3)
                    %if (sum(abs(r(kk,ll,mm,1:3)))~=0 || isnan(sum(r(kk,ll,mm,1:3))) == 0)
                    if (abs(squeeze(r(kk,ll,mm,1))) > 1e-12 || abs(squeeze(r(kk,ll,mm,2))) > 1e-12 || abs(squeeze(r(kk,ll,mm,3))) > 1e-12)
                        % find origin
                        temp_coor=squeeze(r(kk,ll,mm,1:3))';
                        crn_temp=temp_coor-y_unit*half_res-x_unit*half_res-z_unit*half_res;
                        
                        hold on;
                        cube_plot([crn_temp(1),crn_temp(2),crn_temp(3)],res,res,res,'b');
                    end
                end
            end
        end
        
        grid on;xlabel('x');ylabel('y');zlabel('z');
        set(gca,'FontSize',24);set(gca,'FontName','Times New Roman');axis tight; %axis equal;
        view(-45,25)
        
    end
end

