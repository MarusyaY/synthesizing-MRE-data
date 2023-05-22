clear all
close all

%%
tic

disp('Started...');

addpath('./Tools_for_Nifti_and_Analyze_image') %This is to add a matlab library for reading and handling nifti files

%select subject
num_im = num2str(1);

%download the data (no data provided) (example path)
V=load_nii(['./Data/sub' num_im '_50Hz_shear_real.nii.gz']); %real part of complex shear modulus
DR=load_nii(['./Data/sub' num_im '_50Hz_shear_imag.nii.gz']); %imaginary part of complex shear modulus
OR_Re=load_nii('./Data/path'); %reference displacement - real part
OR_Im=load_nii('./Data/path'); %reference displacement - imaginary part

%set according to the data description
dx = 1.5e-3;        % grid point spacing in the x direction [m]
dy = 1.5e-3;        % grid point spacing in the y direction [m]
dz = 1.5e-3;        % grid point spacing in the z direction [m]

%slice to plo
sliceselect=40;

%crop brain
donwsampx1=34; donwsampx2=127;
donwsampy1=20; donwsampy2=153;

Slice=V.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,:);%sliceselect);
Slice_im=DR.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,:);%sliceselect);

orig_re_1=OR_Re.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,:,1);
%orig_im=OR_Im.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,:,1);
orig_re_2=OR_Re.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,:,2);
orig_re_3=OR_Re.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,:,3);


size_of_img = size(V.img);

%size extraction
Nx=size(Slice,1);
Ny=size(Slice,2);
Nz=size(Slice,3);

%set 3D mesh
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);


% %%show selected slice
% figure(1); clf
% contourf(curl_re_1(:,:,sliceselect)); colorbar
% title('orig 1 real');
% 
% figure(2); clf
% contourf(curl_re_2(:,:,sliceselect),20); colorbar
% title('orig 2 real');
% 
% figure(3); clf
% contourf(curl_re_3(:,:,sliceselect)); colorbar
% title('jrig 3 real');


%% brain mask
mask = zeros(Nx,Ny,Nz);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            if Slice(i,j,k)==0
                mask(i,j,k)=0;
            else
                mask(i,j,k)=1;
            end
        end
    end
end

%% constants
v_0 = 1.5; %speed outcide the brain
source_freq = 50; %frequency of the driving force
angular_freq = 2*pi*source_freq;
ro = 1040; %dencity


%inner border of the brain 3D
in_border = zeros(Nx,Ny,Nz);
for i=2:Nx-1
    for j=2:Ny-1
        for k=2:Nz-1
            if (mask(i,j-1,k)==0 || mask(i,j+1,k)==0) && mask(i,j,k)==1
                in_border(i,j,k)=1;
            elseif (mask(i-1,j,k)==0 || mask(i+1,j,k)==0) && mask(i,j,k)==1
                in_border(i,j,k)=1;
            elseif (mask(i,j,k-1)==0 || mask(i,j,k+1)==0) && mask(i,j,k)==1
                in_border(i,j,k)=1;
            else
                in_border(i,j,k)=0;
            end
        end
    end
end

%mask of the external vibration source
F0=in_border(:,:,:);
F0(:,25:end,:)=0;

%% define the properties of the propagation medium

medium.density=ro * ones(Nx, Ny, Nz);

%% attenuation
% 
% atten = zeros(Nx,Ny,Nz);
% for i=1:Nx
%     for j=1:Ny
%         if mask(i,j)==0
%             atten(i,j) = 100;  % [dB/(MHz^y cm)]
%         else
%             atten(i,j) = 0.75;
%         end
%     end
% end
% 
%  medium.alpha_coeff = atten;
%  medium.alpha_power = 1.5;

%% set speed matrix
speed_matrix = zeros(Nx,Ny,Nz);
G_star_mod = sqrt(Slice.^2+Slice_im.^2);
for i = 1:Nx
    for j= 1:Ny
        for k=1:Nz
            if (G_star_mod(i,j,k)==0)
                speed_matrix(i,j,k) = v_0;
            else
                speed_matrix(i,j,k) = sqrt((2*G_star_mod(i,j,k)^2)/((Slice(i,j,k)+G_star_mod(i,j,k))*ro));
            end
        end
    end
end

medium.sound_speed = speed_matrix;  % [m/s]

%% create time array
t_end = 40/60;%2/60;%6e-6; 
dt = 0.0005;% [s]
kgrid.setTime(round(t_end / dt) + 1, dt);

t_start_w = 1250; %for recording
sensor.record_start_index = t_start_w;

%set source
source.p_mask = F0; %border;


source_mag = 1; 
source.p = source_mag * sin(2 * pi * source_freq * kgrid.t_array);

% filter the source to remove any high frequencies not supported by the grid
%source.p = filterTimeSeries(kgrid, medium, source.p);

%% define a time varying sinusoidal source (disk second option)
% karray = kWaveArray;
% 
% p_x = 0; p_y=-67; p_z=0;
% radius = 120;
% diameters = 100;
% 
% f_x =0; f_y=0; f_z=0;
% 
% karray.addBowlElement([p_x, p_y,p_z]*dx, radius*dx, diameters*dx, [f_x,f_y,f_z]*dx);
% 
% source.p_mask = karray.getArrayBinaryMask(kgrid);
% 
% source_mag_1 = 1;     % [Pa]
% sig1 = source_mag_1 * sin(angular_freq * kgrid.t_array);
% source_signal(1, 1:length(sig1)) = sig1;
% 
% source.p = karray.getDistributedSourceSignal(kgrid, source_signal);


%% create a sensor mask covering the entire computational domain using the
% opposing corners of a rectangle
sensor.mask = [1, 1, 1, Nx, Ny, Nz].';

% set the record mode to capture the final wave-field and the statistics at
% each sensor point
sensor.record = {'u'};

% create a display mask to display the transducer
display_mask = source.p_mask;

% assign the input options
input_args = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});



%% calculate displacement from the velocity field

dt = kgrid.dt;

disp_x = zeros(Nx,Ny,Nz,size(sensor_data.ux,4));
disp_y = zeros(Nx,Ny,Nz,size(sensor_data.uy,4));
disp_z = zeros(Nx,Ny,Nz,size(sensor_data.uz,4));


counter = 3;
for i=1:size(disp_x,4)-3
    disp_x(:,:,:,counter)=disp_x(:,:,:,counter-2) + 2*dt*sensor_data.ux(:,:,:,counter-1);
    disp_y(:,:,:,counter)=disp_y(:,:,:,counter-2) + 2*dt*sensor_data.uy(:,:,:,counter-1);
    disp_z(:,:,:,counter)=disp_z(:,:,:,counter-2) + 2*dt*sensor_data.uz(:,:,:,counter-1);
    counter = counter+1;
end

% t_p_start = 700;
% t_p_end = size(sensor_data.ux,4)-3;


%% fit to sine curve
disp_x_final = fit_to_sine_3D(disp_x,angular_freq,0,dt).*mask*1e6;
disp_y_final = fit_to_sine_3D(disp_y,angular_freq,0,dt).*mask*1e6;
disp_z_final = fit_to_sine_3D(disp_z,angular_freq,0,dt).*mask*1e6;
%% curl

[x,y,z] = meshgrid(1:Ny,1:Nx,1:Nz);
Fx = (disp_x_final);
Fy = (disp_y_final);
Fz = (disp_z_final);

[curlx,curly,curlz,cav] = curl(x,y,z,Fx,Fy,Fz);


toc
%%
figure(112);
contourf(real(curlx(:,:,sliceselect)),40); colorbar
title('curlx final');

%%

figure(113);
contourf(real(curly(:,:,sliceselect)),40); colorbar
title('curly final');

%%
figure(114);
contourf(real(curlz(:,:,sliceselect)),20); colorbar
title('curlz final');

%%
figure(115);
contourf(real(disp_x_final(:,:,sliceselect)),20); colorbar
title('x final');

%%

figure(116);
contourf(real(disp_y_final(:,:,sliceselect)),20); colorbar
title('y final');

%%
figure(117);
contourf(real(disp_z_final(:,:,sliceselect)),20); colorbar
title('z final');

%% Evaluation
bias = 0;%10;

%select simulated and reference image
sim = -real(curlz(:,:,:)).*mask*4;
ref = (double(squeeze(orig_re_3(:,:,:))+bias)).*mask;

max_ref = max(max(max(ref)));
min_ref = min(min(min(ref)));

%% plot real and imaginary parts
figure(107); clf
contourf(sim(:,:,sliceselect),30); colorbar
%caxis([min(min(real(sim))),max(max(real(sim)))]);
title('Simulated image');

%%
figure(108); clf
contourf(ref(:,:,sliceselect),10); colorbar
%caxis([min(min(imag(u_final))),max(max(imag(u_final)))]);
title('Reference image');

%% Mean square error

number = 0;
RMSE_re = 0;
for i=1:size(ref,1)
    for j=1:size(ref,2)
        for k = 1:size(ref,3)
            if ref(i,j,k)~=0
                RMSE_re = RMSE_re + (ref(i,j,k)-sim(i,j,k))^2;
                number = number +1;
            end
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
        for k=1:size(ref,3)
            if ref(i,j,k)~=0
                MAE_re = MAE_re + abs(ref(i,j,k)-sim(i,j,k));
                number = number +1;
            end
        end
    end
end

MAE_re = MAE_re/number;
%MAE_pers=100*(MAE_re)/(max_ref-min_ref);

disp(['MAE for real part = ' num2str(MAE_re)]);
%disp(['MAE % = ' num2str(MAE_pers)]);

%% Plot orig vs sim plot along x (real)

figure(10002); clf
hold on
scatter(1:1:Ny, squeeze(sim(30,:,sliceselect)));
scatter(1:1:Ny, squeeze(ref(30,:,sliceselect)));
%ylim([-5 10])
legend('sim','orig')
title('displacement along x-axis for 30th row');
hold off


%% SSIM (there is no toolbox on server for that)
sim = mat2gray(sim);
ref = mat2gray(ref);

figure(1010);
imshow(sim(:,:,45));

figure(1011);
imshow(ref(:,:,45));

function u_final = fit_to_sine_3D(UU,angular_freq,phase0,dt)
%fit_to_sine fits 2D time series target frequency sine curve using fft
% 
%Input: 
%   UU           : 2D time series
%   
%Output:
%   u_final      : matrix, where each pixel is a sine curve 
%                  in complex representation 
Nx = size(UU,1); Ny = size(UU,2); Nz = size(UU,3);
u_final = zeros(Nx,Ny,Nz);



for i=1:Nx
    for j=1:Ny
        for k = 1:Nz
            pointvec=squeeze(UU(i,j,k,:)); 
            N=length(pointvec);
            t = 0:dt:N*dt-dt;
            template = exp(sqrt(-1) * t*angular_freq+phase0);
    
            corr = 2 / N * dot(template,pointvec);
    
            R = abs(corr);
            phi = angle(corr);
    
                    
            u_final(i,j,k) = R*(cos(phi)+sqrt(-1)*sin(phi));
        end
    end
end
end