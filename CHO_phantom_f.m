clear all
close all

tic %records the current time

%% define the size of the object
Nx = 80;
Ny = 80;

disp(['Nx = ' num2str(Nx) ' Ny = ' num2str(Ny)])

%% define the template of the medium distribution inside the object
mask = zeros(Nx,Ny);
counter = 1;
for i=1:Nx
    for j=1:Ny
        if j>=counter
            mask(i,j)=0;
        else
            mask(i,j)=1;
        end
    end
    counter = counter+1;
end

%% set variables
source_freq = 200; %[Hz] frequency of the external driving force
v_1 = 1.9; %[m/s] first velocity
v_2 = 3.5; %[m/s] second velocity

angular_freq = 2*pi*source_freq; %angular frequency
force_amp = 1; %force magnitude
dx = 1.5e-3; %[m] distance between neighboring oscillators

%% Force column vector
%force plate is located in the upper part of the phantom

F0=zeros(Nx,Ny);
source_p=ceil(Ny/3);
F0(Nx,source_p:Ny-source_p)=1;

F0vec=force_amp*reshape(F0,[Nx*Ny,1]);

%% Attenuation

atten=zeros(Nx,Ny)+2;
atten(:,1)=1000;
atten(:,Ny)=1000;
atten(1,:)= 2200 + (4400-2200) .* rand(Ny,1);

% reshape Gamma into diagonal matrix form
atten = reshape(atten,[Nx*Ny,1]);
gamma_attenuation=zeros(Nx*Ny);
for i=1:Nx*Ny
    gamma_attenuation(i,i) = atten(i);
end

%% W - matrix of coupling (force constants)

% calculate k^2 value for each pixel (two mediums)
k_2 = zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        if mask(i,j)==0
            k_2(i,j)= (v_1/dx)^2;
        else
            k_2(i,j)= (v_2/dx)^2;
        end
    end
end

%fill the matrix (see function below)
W = W_matrix_create(k_2);

%% CHO simulation
dt = 1e-4; %time step
t_end = 0.4e-2; %duration of the simulation 

%perform simulation (see function below)
UU = CHO_simulation(angular_freq, Nx, Ny, W, gamma_attenuation, F0vec, dt, t_end);

toc %calculate elapsed time


%% Fit time series to the sine curve (see function below)
u_final = fit_to_sine(UU, angular_freq, 0, dt);
%u_final = fit_to_sine_fft(UU);

%% plot real part
figure(107); clf
contourf(-real(u_final(10:end-20,10:end-10)),20); 
c = colorbar;
c.FontSize = 12;
xlabel('Number of pixels'); 
ylabel('Number of pixels'); 
title('Real part of displacement');

%% plot imaginary part
figure(108); clf
contourf(imag(u_final(10:end-20,10:end-10)),10);
c = colorbar;
c.FontSize = 12;
xlabel('Number of pixels'); 
ylabel('Number of pixels'); 
title('Imaginary part of displacement');

%%
function W = W_matrix_create(k_2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%W_matrix_create creates a coupling matrix for 2D image ready to use in
%  CHO simulation
%
%Input: 
%   k_2 - matrix with k^2 values (force constants)
%   
%Output:
%   W - counping matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx = size(k_2,1);
Ny = size(k_2,2);
W=zeros(Nx*Ny,Nx*Ny); %2D
%filling nth order diagonal (above and under main diagonal)
%coupling between columns
counter = 0;
for j=1:(Ny-1) 
    for i=1:Nx
        %above main diagonal
        W(i+counter,i+Nx+counter)=-k_2(i,j); %[Hz^2] 

        %under main diagonal
        W(i+Nx+counter,i+counter)=-k_2(i,j);
    end
    counter = counter + Nx;
end
%filling 1th order diagonal (above and under main diagonal)
%coupling between rows
counter = 0;
for j=1:Ny
    for i=1:Nx-1
        %above main diagonal
        W(i+counter,i+1+counter)=-k_2(i,j); %[Hz^2] 

        %under main diagonal
        W(i+1+counter,i+counter)=-k_2(i,j);
    end
    counter = counter + Nx;
end
%filling main diagonal
neg_sum_w = -sum(W,2);
for i=1:Nx*Ny
        W(i,i)=neg_sum_w(i);
end

end

function UU = CHO_simulation(angular_freq, Nx, Ny, W, gamma_atten, F0vec, dt, t_end)
%CHO_simulation performes simulation of wave propagation using
%  Coupled Harmonic Oscillator method
%Input: 
%   angular_freq : angular frequency of external driving force
%   Nx           : number of pixels in x direction for original image
%   Ny           : number of pixels in y direction for original image
%   W            : coupling matrix
%   gamma_atten  : attenuation diagonal matrix
%   F0vec        : coulumn vector which denote area and magnitude of
%                  external driving force
%   dt           : time step
%   t_end        : duration of the simulation
%   
%Output:
%   UU           : time series of displacement perpendicular to simulation 
%                  plane 

% Preparation
N =size(W,1);

invAA = inv(W-eye(N)*angular_freq^2);
invBB = inv((W-eye(N)*angular_freq^2) + (angular_freq^2)*gamma_atten*invAA*gamma_atten);
S = angular_freq*invAA*gamma_atten*invBB;

% Main loop

UU = zeros(Nx,Ny,ceil(t_end/dt));
count_image = 1;
phase0=0;

for i=0:dt:t_end
    %solution of differetial equation
    UU_pos_dep = (invBB*F0vec*cos(i*angular_freq+phase0)+S*F0vec*sin(i*angular_freq+phase0))*1e6;
    
    UU(:,:,count_image) = reshape(UU_pos_dep,[Nx,Ny]);

    count_image = count_image + 1;    
end
end


function u_final = fit_to_sine(UU,angular_freq,phase0,dt)
%fit_to_sine fits 2D time series target frequency sine curve using fft
% 
%Input: 
%   UU           : 2D time series
%   angular_freq : angular frequency of external driving force
%   phase0       : phase angle between NMR pulses and sound wave
%   dt           : time step
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