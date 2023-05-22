clear all
close all

tic
% create the computational grid
Nx = 80; % [grid points]
Ny = 80; % [grid points]
dx = 1.5e-3; % [m]
dy = 1.5e-3; % [m]

%set a mesh
kgrid = kWaveGrid(Nx, dx, Ny, dy);

%% define the properties of the propagation medium
v_0 = 20; %[m/s] affects simulation time significantly

v_1 = 1.9; %[m/s]
v_2 = 3.5; %[m/s]
ro = 1040;
source_freq = 200; % [Hz]
angular_freq = 2*pi*source_freq; 
t_end = 7/60; %[s] simulation time


%% select simulation parameters (received from testing)
v_0_options = [10 20 50 80 100 150 200 343 1000];
dt_options = [6e-5 3e-5 1e-5 8e-6 6.5e-6 4.5e-6 3e-06 1.5e-6 6.5e-7];
t_start_w_options = [1000 3000 11000 14000 17000 25000 37000 77000 177000];

[~, ind] = min(abs(v_0_options-v_0));
dt=dt_options(ind);
t_start_w=t_start_w_options(ind);


%% set density
medium.density=ro * ones(Nx, Ny);

%create mask if needed
mask = zeros(Nx,Ny);
mask(1:Nx,1:Ny) = 1;

%%try attenuation
% atten=zeros(Nx,Ny);
% for i=1:Nx
%     atten(i,1:Ny)= 2.2 + (4.4-2.2) .* rand(Ny,1);
% end
% 
% medium.alpha_coeff_shear = atten;  % [dB/(MHz^y cm)] %1.5' %15;
% medium.alpha_coeff_compression = 1.5; %15;
%medium.alpha_power = 1.5;


%create a phantom with two mediums
phantom = ones(Nx,Ny)*v_1;
counter = 1;
for i=1:Nx
    for j=1:Ny
        if j<=counter
            phantom(i,j)=v_2;
        end
    end
    counter = counter+1;
end

phantom = phantom.*mask;

%% set medium properties

medium.sound_speed_shear = phantom;  % [m/s]

speed_comp = v_0*ones(Nx,Ny);

medium.sound_speed_compression = speed_comp;

%% create time array
% t_end = 5/60;%6e-6;       % [s]
% kgrid.makeTime(medium.sound_speed_shear, [], t_end);

kgrid = kWaveGrid(Nx, dx, Ny, dy); % creates the meshgrid 2D
kgrid.setTime(round(t_end / dt) + 1, dt); % create the time array

sensor.record_start_index = t_start_w;


%% define a curved transducer element

F0=zeros(Nx,Ny);
source_p=ceil(Ny/3);
F0(Nx,source_p:Ny-source_p)=1;

source.s_mask = F0;
source_mag_1 = 1;     % [Pa]

source.sxy = source_mag_1 * sin(angular_freq * kgrid.t_array);
source.sxx = 0;
source.syy = 0;

%% create a sensor mask covering the entire computational domain using the
% opposing corners of a rectangle
sensor.mask = [1; 1; Nx; Ny];

% set the record mode to capture the final wave-field and the statistics at
% each sensor point
sensor.record = {'u'};

% create a display mask to display the transducer
display_mask = source.s_mask;

% assign the input options
input_args = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false};

% run the simulation
sensor_data = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});

toc

%%
[x,y] = meshgrid(1:Ny,1:Nx);
dt = kgrid.dt;

disp_x = zeros(Nx,Ny,size(sensor_data.ux,3));
disp_y = zeros(Nx,Ny,size(sensor_data.uy,3));
curl_z = zeros(Nx,Ny,size(sensor_data.ux,3));

counter = 3;
for i=1:size(disp_x,3)-3
    disp_x(:,:,counter)=disp_x(:,:,counter-2) + 2*dt*sensor_data.ux(:,:,counter-1);
    disp_y(:,:,counter)=disp_y(:,:,counter-2) + 2*dt*sensor_data.uy(:,:,counter-1);
    [curlz,cav] = curl(x,y,sensor_data.ux(:,:,counter),sensor_data.uy(:,:,counter));

    curl_z(:,:,counter) = curlz;
    counter = counter+1;
end


%% fit to sine curve

phase0=0;
disp_x_final = fit_to_sine(disp_x,angular_freq,phase0,dt)*1e6;
disp_y_final = fit_to_sine(disp_y,angular_freq,phase0,dt)*1e6;
curl_z_final = fit_to_sine(curl_z,angular_freq,phase0,dt)*1e4;


disp_x_re = real(disp_x_final);%.*mask;
disp_y_re = real(disp_y_final);%.*mask;
curl_z_re = real(curl_z_final);


%% plot real and imaginary parts curl
figure(107); clf
contourf(disp_x_re(10:end-20,10:end-10),10); 
c = colorbar;
c.FontSize = 12;
%caxis([min(min(disp_x_re))*0.1,0.1*max(max(disp_x_re))]);
title('final displacement - real x');

%%
figure(108); clf
contourf(disp_y_re(10:end-20,10:end-10)); 
c = colorbar;
c.FontSize = 12;
%caxis([min(min(imag(u_final))),max(max(imag(u_final)))]);
title('final displacement - real y');


%% curl of real part of displacement
[curlz,cav] = curl(x,y,disp_x_re,disp_y_re);
%%
figure(109); clf
contourf(-curlz(10:end-20,10:end-10)); colorbar;
%caxis([min(min(imag(u_final))),max(max(imag(u_final)))]);
title('curl of final displacement');


disp(['compressional wave velocity = ' num2str(v_0)]);
disp(['t_end = ' num2str(t_end) '; dt = ' num2str(dt) '; t_start_w = ' num2str(t_start_w)]); 
%disp(['Force: [' num2str(f_x_1) ',' num2str(f_y_1) '],[' num2str(f_x_2) ',' num2str(f_y_2) ']; ' tensor]);



function u_final = fit_to_sine(UU,angular_freq,phase0,dt)
%fit_to_sine fits 2D time series target frequency sine curve using fft
% 
%Input: 
%   UU           : 2D time series
%   angular_frequency
%   phase0
%   dt
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

% for i = 1:Nx
%     for j = 1:Ny
%         vect = squeeze(UU(i,j,:));
%         N = length(vect);
%         d_fft=fft(vect); %Fast Fourier Transform
%         [~,I] = max(abs(d_fft(1:ceil(N/2)+1))); %index of max frequency 
%      
%         f = fs*(0:(N/2))/N;
%         f(I)
%         
%         R = 2*abs(d_fft(I)/N); %magnitude 
%         phi = angle(d_fft(I)); %phase
%         
%         % complex representation
%         u_final(i,j) = R*(cos(phi)+sqrt(-1)*sin(phi));
% 
%     end
end
end