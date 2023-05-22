clear all
close all

%%
disp('Started...');
tic

addpath('./Tools_for_Nifti_and_Analyze_image') %This is to add a matlab library for reading and handling nifti files
%addpath('../k-Wave') %could be installed via apps

%chose subject
num_im = num2str(1);

%% constants
source_freq = 50; %[Hz]
angular_freq = 2*pi*source_freq;
ro = 1040;
v_0 = 1.5; %velocity outside the brain

%select plane
plane = 1; % 1 - horizontal plane; 2 - sagittal plane

%select excitation scheme
scheme = 1; % 1- back of the head or 2 - ark source outside the brain

%download the data (no data provided) (example path)
V=load_nii(['./Data/sub' num_im '_50Hz_shear_real.nii.gz']); %real part of complex shear modulus
DR=load_nii(['./Data/sub' num_im '_50Hz_shear_imag.nii.gz']); %imaginary part of complex shear modulus
OR_Re=load_nii('./Data/path'); %reference displacement - real part
OR_Im=load_nii('./Data/path'); %reference displacement - imaginary part

%set according to the data description
dx = 1.5e-3;        % grid point spacing in the x direction [m]
dy = 1.5e-3;        % grid point spacing in the y direction [m]


%plane selection
if plane == 1
    %% xy (horizontal plane)
    sliceselect = 45;
    donwsampx1=34; donwsampx2=127;
    donwsampy1=20; donwsampy2=153;
    
    Slice=V.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect);
    Slice_im=DR.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect);
    orig_re_3=OR_Re.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,3);
    orig_re_2=OR_Re.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,2);
    orig_re_1=OR_Re.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,1);
    orig_im_3=OR_Im.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,3);
    
else
    % xz (sagittal plane)
    sliceselect=80;
    donwsampy1=20; donwsampy2=153;
    Slice=squeeze(V.img(sliceselect,donwsampy1:donwsampy2,:))';
    Slice_im=squeeze(DR.img(sliceselect,donwsampy1:donwsampy2,:))';
  
    orig_re_1=squeeze(OR_Re.img(sliceselect,donwsampy1:donwsampy2,:,1))';
    orig_re_2=squeeze(OR_Re.img(sliceselect,donwsampy1:donwsampy2,:,2))';
    orig_re_3=squeeze(OR_Re.img(sliceselect,donwsampy1:donwsampy2,:,3))';
end

size_of_img = size(V.img);

%size extraction
Nx=size(Slice,1);
Ny=size(Slice,2);

%set mesh
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% %% show reference images
% figure(1); clf
% contourf(orig_re_1); colorbar
% title('orig real 1');
% 
% figure(2); clf
% contourf(orig_re_2); colorbar
% title('orig real 2');
% 
% figure(3); clf
% contourf(orig_re_3); colorbar
% title('orig real 3');


%% brain mask
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


%% define the properties of the propagation medium

medium.density=ro * ones(Nx, Ny);

speed_matrix = zeros(Nx,Ny);
G_star_mod = sqrt(Slice.^2+Slice_im.^2);
for i = 1:Nx
    for j= 1:Ny
        if (mask(i,j)==0)
            speed_matrix(i,j) = v_0; %speed outside the brain
        else
            speed_matrix(i,j) = sqrt((2*G_star_mod(i,j)^2)/((Slice(i,j)+G_star_mod(i,j))*ro));
        end
    end
end

medium.sound_speed = speed_matrix;  % [m/s]

%% create time array
t_end = 20/60;%2/60;%6e-6; 
dt = 0.0002; % [s]
kgrid.setTime(round(t_end / dt) + 1, dt);

% t_end = 50/60;%6e-6;       % [s]
%kgrid.makeTime(medium.sound_speed, [], t_end);

 t_start_w = 900; %for recording
 sensor.record_start_index = t_start_w;


%% define the inner boundary (see function below)

in_border = define_in_border(mask);

% select excitation scheme
if scheme == 1
    %% define force along inner boundary (first option)
    F0=in_border(:,:);
    F0(:,25:end)=0;
    
    source.p_mask = F0;
    source_mag_1 = 1;
    source.p = source_mag_1 * sin(angular_freq * kgrid.t_array);
else
    %% define force (second option)
    karray = kWaveArray;
    
    p_x = 0; p_y=-67;
    radius = 70;
    diameters = 100;
    
    f_x =0; f_y=0;
    
    karray.addArcElement([p_x, p_y]*dx, radius*dx, diameters*dx, [f_x,f_y]*dx);
    
    source.p_mask = karray.getArrayBinaryMask(kgrid);
    
    source_mag_1 = 1;     % [Pa]
    sig1 = source_mag_1 * sin(angular_freq * kgrid.t_array);
    source_signal(1, 1:length(sig1)) = sig1;
    
    source.p = karray.getDistributedSourceSignal(kgrid, source_signal);
end

%% create a sensor mask covering the entire computational domain using the
% opposing corners of a rectangle
sensor.mask = [1; 1; Nx; Ny];

% set the record mode to capture the final wave-field and the statistics at
% each sensor point
sensor.record = {'u'};

% create a display mask to display the transducer
display_mask = source.p_mask;

% assign the input options (here gpu could be added)
input_args = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

toc

%% calculate displacement from estimated velocity field
[x,y] = meshgrid(1:Ny,1:Nx);
dt = kgrid.dt;

disp_x = zeros(Nx,Ny,size(sensor_data.ux,3));
disp_y = zeros(Nx,Ny,size(sensor_data.uy,3));
curl_d = zeros(Nx,Ny,size(sensor_data.ux-3,3));

counter = 1;
for i=1:size(disp_x,3)-3
    if i>2    
        disp_x(:,:,counter) = disp_x(:,:,counter-2) + 2*dt*sensor_data.ux(:,:,counter-1);
        disp_y(:,:,counter) = disp_y(:,:,counter-2) + 2*dt*sensor_data.uy(:,:,counter-1);
    end
    [curl_d(:,:,counter),cav]=curl(x,y,disp_x(:,:,counter),disp_y(:,:,counter));
    counter = counter+1;
end
%%




%% gif
% t_p_start = 1;
% t_p_end = size(disp_x,3)-3;
% counter = 1;
% for i=t_p_start:size(disp_x,3)-3
%     fig = figure(200); clf
%     contourf(real(disp_x(:,:,i)).*mask); colorbar
%     title('wave propagation');
%     minval=min(min(disp_x(:,:,i)))
%     caxis([minval, max(max(disp_x(:,:,i)))])
%     pause(0.005)
%     drawnow
%     frame = getframe(fig);
%     im{counter} = frame2im(frame);
%     counter = counter + 1;
% end

% %% to save gif file
% 
% filename = "./Save/wave_propagation.gif"; % Specify the output file name
% for idx = 1:size(im,2)
%     [L,map] = rgb2ind(im{idx},256);
%     if idx == 1
%         imwrite(L,map,filename,"gif","LoopCount",Inf,"DelayTime",2);
%     else
%         imwrite(L,map,filename,"gif","WriteMode","append","DelayTime",2);
%     end
% end

%%
figure(333); clf
hold on
%new=reshape(v(50,50,:),3334);
for  i=1:size(curl_d,3)-1
    plot(i,curl_d(45,64,i),'x');
end
hold off
%%
figure(444);clf
hold on
for i=1:20:Nx
    plot(curl_d(i,:,end-4));
end
hold off

%% Fit to the sinusoid

disp_x_final = fit_to_sine(disp_x(:,:,1:end-2),angular_freq,0,dt);
disp_y_final = fit_to_sine(disp_y(:,:,1:end-2),angular_freq,0,dt);
curl_d_final = fit_to_sine(curl_d(:,:,1:end-2),angular_freq,0,dt);


disp_x_final=disp_x_final.*mask*1e6;
disp_y_final=disp_y_final.*mask*1e6;
curl_d_final = curl_d_final.*mask*1e6*4;


%% Evaluation
bias = 0;

%set simulation and reference images
if plane == 1
    ref = double(orig_re_3+bias).*mask;
    sim = -real(curl_d_final); 
else
    ref = double(orig_re_1+bias).*mask;
    sim = -real(curl_d_final);
end

max_ref = max(max(ref));
min_ref = min(min(ref));

%% plot real and imaginary parts
figure(107); clf
contourf(sim,20); colorbar
%caxis([min(min(real(u_final))),max(max(real(u_final)))]);
title('Simulated image');

%%
figure(108); clf
contourf(ref,15); colorbar
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
ylim([-6 6])
lgd = legend('Data','Fit');

%% SSIM (there is no toolbox on server for that)
sim = mat2gray(sim);
ref = mat2gray(ref);

figure(109);
imshow(sim);

figure(110);
imshow(ref);

SSIM_re = ssim(sim,ref);
disp(['SSIM for real part = ' num2str(SSIM_re)]);


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

function u_final = fit_to_sine_fft(UU)
%fit_to_sine fits 2D time series to sime curve using fft
% 
%Input: 
%   UU           : 2D time series
%   
%Output:
%   u_final      : matrix, where each pixel is a sine curve 
%                  in complex representation 

Nx = size(UU,1); Ny = size(UU,2);
u_final = zeros(Nx,Ny);

for i = 1:Nx
    for j = 1:Ny
        vect = squeeze(UU(i,j,:));
        N = length(vect);
        d_fft=fft(vect); %Fast Fourier Transform
        [~,I] = max(abs(d_fft(1:ceil(N/2)+1))); %index of max frequency 
        
        R = 2*abs(d_fft(I)/N); %magnitude 
        phi = angle(d_fft(I)); %phase
        
        % complex representation
        u_final(i,j) = R*(cos(phi)+sqrt(-1)*sin(phi));
    end
end
end