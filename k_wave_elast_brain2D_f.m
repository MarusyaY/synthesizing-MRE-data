clear all
close all

disp('Started...');

addpath('./Tools_for_Nifti_and_Analyze_image') %This is to add a matlab library for reading and handling nifti files

%chose subject
num_im = num2str(1);

%select simulation plane
plane = 1; % 1 - horizontal plane; 2 - sagittal plane

%select excitation scheme
scheme = 1; % 1 - back of the head or 2 - arc source outside the brain

%select stress tensor combibation
comb = 3; % 1 sxx, 2 syy, 3 sxy, 4 sxx&syy 5 sxx&sxy 6 syy&sxy, 7 all 

%% constants
v_comp = 20; %velocity of compressional waves
v_mean_brain = 1.5; % velocity of shear wave outside the brain
ro = 1040; %density of tissues
source_freq = 50; % [Hz] frequency of the external driving force
angular_freq = 2*pi*source_freq;
t_end = 7/60; %simulation time

%% set simulation parameters (received from testing)
v_0_options = [10 20 50 80 100 150 200 343 1000];
dt_options = [6e-5 3e-5 1e-5 8e-6 6.5e-6 4.5e-6 3e-06 1.5e-6 6.5e-7];
t_start_w_options = [1000 3000 11000 14000 17000 25000 37000 77000 177000];

[~, ind] = min(abs(v_0_options-v_comp));
dt=dt_options(ind); %time step
t_start_w=t_start_w_options(ind); %when data recording begins

%download the data (no data provided) (example path)
V=load_nii(['./Data/sub' num_im '_50Hz_shear_real.nii.gz']); %real part of complex shear modulus
DR=load_nii(['./Data/sub' num_im '_50Hz_shear_imag.nii.gz']); %imaginary part of complex shear modulus
OR_Re=load_nii('./Data/path'); %reference displacement - real part
OR_Im=load_nii('./Data/path'); %reference displacement - imaginary part

%set according to the data description
dx = 1.5e-3;        % grid point spacing in the x direction [m]
dy = 1.5e-3;        % grid point spacing in the y direction [m]

%set simulation plane
if plane == 1
    % xy (horizontal plane)
    sliceselect=40;
    donwsampx1=34; donwsampx2=127;
    donwsampy1=20; donwsampy2=153;
    
    Slice=V.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect);
    Slice_im=DR.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect);
    orig_re_3=OR_Re.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,3);
    orig_re_2=OR_Re.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,2);
    orig_re_1=OR_Re.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,1);
    orig_im_3=OR_Im.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,3);
    orig_im_2=OR_Im.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,2);
    orig_im_1=OR_Im.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,1);
    
else
    %% xz (sagittal plane)
    sliceselect=80;
    donwsampy1=20; donwsampy2=153;
    Slice=squeeze(V.img(sliceselect,donwsampy1:donwsampy2,:))';
    Slice_im=squeeze(DR.img(sliceselect,donwsampy1:donwsampy2,:))';
    
    orig_re_1=squeeze(OR_Re.img(sliceselect,donwsampy1:donwsampy2,:,1))';
    %curl_im=Curl_im.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,1);
    orig_re_2=squeeze(OR_Re.img(sliceselect,donwsampy1:donwsampy2,:,2))';
    orig_re_3=squeeze(OR_Re.img(sliceselect,donwsampy1:donwsampy2,:,3))';
end

%% define the size of the image
Nx = size(Slice,1);
Ny = size(Slice,2)


%% mask of the brain

mask = zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        if Slice(i,j)==0
            mask(i,j)=0;
        else
            mask(i,j)=1;
        end
    end
end

%% time array

kgrid = kWaveGrid(Nx, dx, Ny, dy); % creates the meshgrid 2D
kgrid.setTime(round(t_end / dt) + 1, dt); % create the time array

sensor.record_start_index = t_start_w;

%% Define the properties of the propagation medium
mask_special = mask;
mask_special(1:end,1:end)=1;

medium.density = ro * ones(Nx, Ny);
medium.sound_speed_compression = ones(Nx,Ny)*v_comp;


speed_matrix = zeros(Nx,Ny);
G_star_mod = sqrt(Slice.^2+Slice_im.^2);
for i = 1:Nx
    for j= 1:Ny
        if (mask(i,j)~=0)
            speed_matrix(i,j) = sqrt((2*G_star_mod(i,j)^2)./((Slice(i,j)+G_star_mod(i,j))*ro));
        else 
            speed_matrix(i,j)=v_mean_brain;
        end
    end
end

medium.sound_speed_shear = speed_matrix; 
 
%% define the inner boundary (just in case, not using now)

in_border = define_in_border(mask);

%set exitation scheme
if scheme == 1
    % define force along inner boundary (first option)
    
    F0=in_border(:,:);
    F0(:,25:end)=0;
    
    source.s_mask = F0;
    source_mag_1 = 1; %0.5; % [Pa]
    
    source_signal = source_mag_1 * sin(angular_freq * kgrid.t_array);
    source_signal_0 = 0;
    
    %set stress tensor (see function below)
    [sourse.sxx, source.syy, source.sxy] = define_tensor(comb, source_signal, source_signal_0);
    
    tensor = 'sxy';
else
    %% define force (second option)
    karray = kWaveArray;
    
    p_x = 0; p_y=-67;
    radius = 70;
    diameters = 70;
    
    f_x =0; f_y=0;
    
    karray.addArcElement([p_x, p_y]*dx, radius*dx, diameters*dx, [f_x,f_y]*dx);
    %create arc source
    source.s_mask = karray.getArrayBinaryMask(kgrid);
    
    source_mag_1 = 1; 
    sig1 = source_mag_1 * sin(angular_freq * kgrid.t_array);
    source_signal(1, 1:length(sig1)) = sig1;
    
    source_s = karray.getDistributedSourceSignal(kgrid, source_signal);
    source_s_0 = karray.getDistributedSourceSignal(kgrid, 0);
     %set stress tensor (see function below)
    [sourse.sxx, source.syy, source.sxy] = define_tensor(comb, source_s, source_s_0);

end
% where to record
sensor.mask = [1; 1; Nx; Ny];

sensor.record ={'u'}; %Sets what to record at each voxel

display_mask=source.s_mask;

% assign the input options
input_args = {'DisplayMask', display_mask, 'DataCast','single', 'PMLInside', false, 'PlotPML', false, 'PMLAlpha',[4,4]};

%simulation
sensor_data = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});


%% post-processinng

%calculate displacement from velocity field
[x,y] = meshgrid(1:Ny,1:Nx);
dt = kgrid.dt;

disp_x = zeros(Nx,Ny,size(sensor_data.ux,3));
disp_y = zeros(Nx,Ny,size(sensor_data.uy,3));

counter = 3;
for i=1:size(disp_x,3)-3
    disp_x(:,:,counter)=disp_x(:,:,counter-2) + 2*dt*sensor_data.ux(:,:,counter-1);
    disp_y(:,:,counter)=disp_y(:,:,counter-2) + 2*dt*sensor_data.uy(:,:,counter-1);

    counter = counter+1;
end

%% fit to sine curve
disp_x_final = fit_to_sine(disp_x,angular_freq,0,dt).*mask*1e6;


disp_y_final = fit_to_sine(disp_y,angular_freq,0,dt).*mask*1e6;


%% input values
disp(['compressional wave velocity = ' num2str(v_comp)]);
disp(['t_end = ' num2str(t_end) '; dt = ' num2str(dt) '; t_start_w = ' num2str(t_start_w)]); 


%% Evaluation
bias = 10; %7/8/0; 2sub%4.5/8; 1sub7/10

%select simulation and reference data 
sim = -3*double(real(disp_y_final))/10;
ref = (double(orig_re_2)+bias).*mask;

max_ref = max(max(ref));
min_ref = min(min(ref));

%% plot real and imaginary parts
figure(107); clf
contourf(sim,20); colorbar
%caxis([min(min(real(u_final))),max(max(real(u_final)))]);
title('Simulated image');

%%
figure(108); clf
contourf(ref,10); colorbar
%caxis([min(min(imag(u_final))),max(max(imag(u_final)))]);
title('Reference image');

%% Mean square error

number = 0;
RMSE_re = 0;
for i=1:size(ref,1)
    for j=1:size(ref,2)
        if ref(i,j)~=0
            RMSE_re = RMSE_re + (ref(i,j)-sim(i,j))^2;
            number = number +1;
        end
    end
end

RMSE_re = sqrt(RMSE_re/number);

disp(['RMSE for real part = ' num2str(RMSE_re) '; Max ref value = ' num2str(max_ref) '; Min sim value = ' num2str(min_ref)]);

%% mean absolute error

number = 0;
MAE_re = 0;
for i=1:size(ref,1)
    for j=1:size(ref,2)
        if ref(i,j)~=0
            MAE_re = MAE_re + abs(ref(i,j)-sim(i,j));
            number = number +1;
        end
    end
end

MAE_re = MAE_re/number;
MAE_pers=100*(MAE_re)/(max_ref-min_ref);


disp(['MAE for real part = ' num2str(MAE_re) ' MAE % = ' num2str(MAE_pers)]);


%% Plot orig vs sim plot along x (real)

figure(10002); clf
hold on
scatter(1:1:Ny, squeeze(sim(50,:)));
scatter(1:1:Ny, squeeze(ref(50,:)));
%ylim([-5 10])
legend('sim','orig')
title('displacement along x-axis for 50th row');
hold off

%% plot simulation vs reference
sim_reshape = reshape(sim,[Nx*Ny,1]);
ref_reshape = reshape(ref,[Nx*Ny,1]);

b1=ref_reshape\sim_reshape;
lin = ref_reshape*b1;
X = [ones(length(ref_reshape),1) ref_reshape];
b2 = X\sim_reshape;
lin2 = X*b2;

figure(1010);
clf
plot (ref_reshape, sim_reshape,'o');
hold on
%plot (ref_reshape, lin);
plot (ref_reshape, lin2);
xlabel('Reference displacement', 'fontsize',14);
ylabel('Simulated displacement', 'fontsize',14);
xlim([-7 5])
lgd = legend('Data','Fit');

%% SSIM (there is no toolbox on server for that)
sim = mat2gray(sim);
ref = mat2gray(ref);

figure(1109);
imshow(sim);

figure(1110);
imshow(ref);

SSIM_re = ssim(sim,ref);
disp(['SSIM for real part = ' num2str(SSIM_re) ' bias = ' num2str(bias)]);


function in_border = define_in_border(mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%in_border create a mask of inner boundary of the brain
%
%Input: 
%   mask - mask of the brain
%   
%Output:
%   in_border - mask of the inner boundary of the brain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = size(mask,1); Ny = size(mask,2);
in_border = zeros(Nx,Ny);
for i=2:Nx-1
    for j=2:Ny-1
        if (mask(i,j-1)==0 || mask(i,j+1)==0) && mask(i,j)==1
            in_border(i,j)=1;
        elseif  (mask(i-1,j)==0 || mask(i+1,j)==0) && mask(i,j)==1
            in_border(i,j)=1;
        else
            in_border(i,j)=0;
        end
    end
end
end


function u_final = fit_to_sine(UU,angular_freq,phase0,dt)
%fit_to_sine fits 2D time series target frequency sine curve using fft
% 
%Input: 
%   UU           : 2D time series
%   
%Output:
%   u_final      : matrix, where each pixel is a sine curve 
%                  in complex representation 
Nx = size(UU,1); Ny = size(UU,2);
u_final = zeros(Nx,Ny);



for i=1:Nx
    for j=1:Ny               
        pointvec=squeeze(UU(i,j,:)); 
        N=length(pointvec);
        t = 0:dt:N*dt-dt;
        template = exp(sqrt(-1) * t*angular_freq+phase0);

        corr = 2 / N * dot(template,pointvec);

        R = abs(corr);
        phi = angle(corr);

                
        u_final(i,j) = R*(cos(phi)+sqrt(-1)*sin(phi));
          
    end
end
end


function [sxx, syy, sxy] = define_tensor(comb, source_signal, source_signal_0)
if comb == 1
    sxx = source_signal;
    syy = source_signal_0;
    sxy = source_signal_0;
elseif comb == 2
    sxx = source_signal_0;
    syy = source_signal;
    sxy = source_signal_0;
elseif comb == 3
    sxx = source_signal_0;
    syy = source_signal_0;
    sxy = source_signal;
elseif comb == 4
    sxx = source_signal;
    syy = source_signal;
    sxy = source_signal_0;
elseif comb == 5
    sxx = source_signal;
    syy = source_signal_0;
    sxy = source_signal;
elseif comb == 6
    sxx = source_signal_0;
    syy = source_signal;
    sxy = source_signal;
else
    sxx = source_signal;
    syy = source_signal;
    sxy = source_signal;
end
end