clear all
close all

tic
% create the computational grid
Nx = 80; % [grid points]
Ny = 80; % [grid points]
dx = 1.5e-3; % [m]
dy = 1.5e-3; % [m]

%create mesh
kgrid = kWaveGrid(Nx, dx, Ny, dy);

%% define the properties of the propagation medium
v_1 = 1.9; %[m/s]
v_2 = 3.5; %[m/s]
ro = 1040;

medium.density=ro * ones(Nx, Ny);

phantom = zeros(Nx,Ny);
counter = 1;
for i=1:Nx
    for j=1:Ny
        if j<=counter
            phantom(i,j)=v_2;
        else
            phantom(i,j)=v_1;
        end
    end
    counter = counter+1;
end

medium.sound_speed = phantom;  % [m/s]


%% create time array
t_end = 6/60;%6e-6;       % [s]
kgrid.makeTime(medium.sound_speed, [], t_end);

%% define a transducer element

F0=zeros(Nx,Ny);
source_p=ceil(Ny/3);
F0(Nx,source_p:Ny-source_p)=1;

source.p_mask = F0;

%% define a time varying sinusoidal source
source_freq = 200; %0.25e6;       % [Hz]
source_mag = 1; %0.5;           % [Pa]
source.p = source_mag * sin(2 * pi * source_freq * kgrid.t_array);

% filter the source to remove any high frequencies not supported by the grid
source.p = filterTimeSeries(kgrid, medium, source.p);

%% create a sensor mask covering the entire computational domain using the
% opposing corners of a rectangle
sensor.mask = [1; 1; Nx; Ny];

% set the record mode to capture the final wave-field and the statistics at
% each sensor point
sensor.record = {'u'};

% create a display mask to display the transducer
display_mask = source.p_mask;

% assign the input options
input_args = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

toc

%% calculate displacement from velocities
[x,y] = meshgrid(1:Ny,1:Nx);
dt = kgrid.dt;

disp_x = zeros(Nx,Ny,size(sensor_data.ux,3));
disp_y = zeros(Nx,Ny,size(sensor_data.uy,3));
curl_z = zeros(Nx,Ny,size(sensor_data.ux-3,3));

counter = 3;
for i=1:size(disp_x,3)-3
    disp_x(:,:,counter)=disp_x(:,:,counter-2) + 2*dt*sensor_data.ux(:,:,counter-1);
    disp_y(:,:,counter)=disp_y(:,:,counter-2) + 2*dt*sensor_data.uy(:,:,counter-1);
    [curlz,cav] = curl(x,y,disp_x(:,:,counter),disp_y(:,:,counter));
    curl_z(:,:,counter) = curlz;
    counter = counter+1;
end




% %% gif
% t_p_start = 1000;
% t_p_end = size(curl_z,3)-3;
% counter = 1;
% for i=t_p_start:t_p_end
%     fig = figure(200); clf
%     contourf(real(curl_z(:,:,i))); colorbar
%     title('wave propagation');
%     minval=min(min(curl_z(:,:,i)))
%     caxis([minval, max(max(curl_z(20:100,20:100,i)))])
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

%% Fit to the sinusoid

curl_final = fit_to_sine_freq(curl_z(:,:,600:end),2*pi*source_freq,0,dt)*1e6;

%% plot real and imaginary parts curl
figure(107); clf
contourf(-real(curl_final(10:end-20,10:end-10)),10);
c = colorbar;
c.FontSize = 12;
colorbar
%caxis([min(min(real(curl_final))),max(max(real(curl_final)))]);
title('k-Wave - kspaceFirstOrder');

%%
% figure(108); clf
% contourf(imag(curl_final(10:end-20,10:end-10))); 
% c = colorbar;
% c.FontSize = 12;
% %caxis([min(min(imag(u_final))),max(max(imag(u_final)))]);
% %title('curl.displacement imaginary');
%%

grey = mat2gray(real(curl_final));
%greq = normalize(grey);
figure(109);
imshow(grey(10:end-20,10:end-10));



function u_final = fit_to_sine_fft(UU,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fit_to_sine_fft fits 2D time series to sime curve using fft
% 
%Input: 
%   UU           : 2D time series
%   
%Output:
%   u_final      : matrix, where each pixel is a sine curve 
%                  in complex representation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = size(UU,1); Ny = size(UU,2);
u_final = zeros(Nx,Ny);

fs=1/dt;

for i = 1:Nx
    for j = 1:Ny
        vect = squeeze(UU(i,j,:));
        N = length(vect);
        d_fft=fft(vect); %Fast Fourier Transform
        [~,I] = max(abs(d_fft(1:ceil(N/2)+1))); %index of max frequency 
        
        f = fs*(0:(N/2))/N; % Frequency vector
        f(I);
        
        R = 2*abs(d_fft(I)/N); %magnitude 
        phi = angle(d_fft(I)); %phase
        
        % complex representation
        u_final(i,j) = R*(cos(phi)+sqrt(-1)*sin(phi));
%         t = 0:dt:N*dt-dt;
%           figure(20000); clf
%           plot(t,cos(2*pi*f(I)*t+phi)*R, 'x')
%           hold on
%           plot(t,vect)
%           pause(1)
%           hold off
    end
end
end

function u_final = fit_to_sine_freq(UU,angular_freq,phase0,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fit_to_sine fits 2D time series target frequency sine curve using fft
% 
%Input: 
%   UU           : 2D time series
%   
%Output:
%   u_final      : matrix, where each pixel is a sine curve 
%                  in complex representation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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