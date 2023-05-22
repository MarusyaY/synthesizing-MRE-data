clear all
close all

%%
tic

disp('Started...');
%This is to add a matlab library for reading and handling nifti files
addpath('./Tools_for_Nifti_and_Analyze_image') 

%chose subject
num_im = num2str(1);

%% variables
source_freq = 50; %[Hz]
angular_freq = 2*pi*source_freq;
force_amp = 1; 
k_0 = 0.01;
ro = 1040;

%simulation plane
% 1 - horizontal plane; 2 - sagittal plane
plane = 1; 

%excitation scheme
% 1 - back of the head 
% 2 - back of the head and front of the head, bur in the opposit direction
scheme = 1; 


%download the data
V=load_nii(['./Data/U01_UDEL_000' num_im '_01_v3/U01_UDEL_000' num_im '_01_MRE_AP_50Hz/U01_UDEL_000' num_im '_01_MRE_AP_50Hz_props_shear_real.nii.gz']);
DR=load_nii(['./Data/U01_UDEL_000' num_im '_01_v3/U01_UDEL_000' num_im '_01_MRE_AP_50Hz/U01_UDEL_000' num_im '_01_MRE_AP_50Hz_props_shear_imag.nii.gz']);
Curl_re = load_nii(['./Data/U01_UDEL_000' num_im '_01_v3/U01_UDEL_000' num_im '_01_MRE_AP_50Hz/U01_UDEL_000' num_im '_01_MRE_AP_50Hz_curl_re.nii.gz']);
Curl_im = load_nii(['./Data/U01_UDEL_000' num_im '_01_v3/U01_UDEL_000' num_im '_01_MRE_AP_50Hz/U01_UDEL_000' num_im '_01_MRE_AP_50Hz_curl_im.nii.gz']);
OR_Re=load_nii(['./Data/U01_UDEL_000' num_im '_01_v3/U01_UDEL_000' num_im '_01_MRE_AP_50Hz/U01_UDEL_000' num_im '_01_MRE_AP_50Hz_disp_re.nii.gz']);
OR_Im=load_nii(['./Data/U01_UDEL_000' num_im '_01_v3/U01_UDEL_000' num_im '_01_MRE_AP_50Hz/U01_UDEL_000' num_im '_01_MRE_AP_50Hz_disp_im.nii.gz']);

%set according to the data description
dx = 1.5e-3;        % grid point spacing in the x direction [m]
dy = 1.5e-3;        % grid point spacing in the y direction [m]

% cropped 2D slice (needs to be adjasted with respect to the dataset)
if plane == 1
%% xy horizontal plane
    sliceselect = 45;
    %chose patch size
    donwsampx1=34; donwsampx2=127;
    donwsampy1=20; donwsampy2=153;
    % cropp images
    Slice=V.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect);
    Slice_im=DR.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect);
    orig_re_3=OR_Re.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,3);
    orig_re_2=OR_Re.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,2);
    orig_re_1=OR_Re.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,1);
    orig_im_3=OR_Im.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,3);
    orig_im_2=OR_Im.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,2);
    orig_im_1=OR_Im.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,1);
    
    curl_re_1=Curl_re.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,1);
    %curl_im=Curl_im.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,1);
    curl_re_2=Curl_re.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,2);
    curl_re_3=Curl_re.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,3);
else
    %% xz (perpendicular to left-right plane)
    sliceselect=80;
    %chose patch size
    donwsampy1=20; donwsampy2=153;
    % cropp images
    Slice=squeeze(V.img(sliceselect,donwsampy1:donwsampy2,:))';
    Slice_im=squeeze(DR.img(sliceselect,donwsampy1:donwsampy2,:))';
     
    curl_re_1=squeeze(Curl_re.img(sliceselect,donwsampy1:donwsampy2,:,1))';
    %curl_im=Curl_im.img(donwsampx1:donwsampx2,donwsampy1:donwsampy2,sliceselect,1);
    curl_re_2=squeeze(Curl_re.img(sliceselect,donwsampy1:donwsampy2,:,2))';
    curl_re_3=squeeze(Curl_re.img(sliceselect,donwsampy1:donwsampy2,:,3))';
    
    orig_re_1=squeeze(OR_Re.img(sliceselect,donwsampy1:donwsampy2,:,1))';
    orig_re_2=squeeze(OR_Re.img(sliceselect,donwsampy1:donwsampy2,:,2))';
    orig_re_3=squeeze(OR_Re.img(sliceselect,donwsampy1:donwsampy2,:,3))';
    
    orig_im_1=squeeze(OR_Im.img(sliceselect,donwsampy1:donwsampy2,:,1))';
    orig_im_2=squeeze(OR_Im.img(sliceselect,donwsampy1:donwsampy2,:,2))';
    orig_im_3=squeeze(OR_Im.img(sliceselect,donwsampy1:donwsampy2,:,3))';
end

size_of_img = size(V.img);

%size extraction
Nx=size(Slice,1);
Ny=size(Slice,2);

disp(['Nx = ' num2str(Nx) ' Ny = ' num2str(Ny)])


%% show reference images
% figure(1); clf
% contourf(orig_re_1); colorbar
% title('orig real 1');
% 
% figure(2); clf
% contourf(orig_re_2); colorbar
% title('orig real 2');
% 
% figure(3); clf
% contourf(orig_re_3);
% cbar = colorbar;
% l = ylabel(cbar,'Displacement (μm)','Rotation',270);
% l.Position(1) =3;
% %caxis([min(min(real(u_final))),max(max(real(u_final)))]);
% xlabel('Number of pixels'); 
% ylabel('Number of pixels'); 
% title('Reference displacement along z axis (top to bottom)');

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

%% define and plot force vector (see function below)
%define inner boundary of the brain
in_border = define_in_border(mask);

%this should be adjasted according to data set parameters
F0 = set_force(scheme,in_border);

F0vec=force_amp*reshape(F0,[Nx*Ny,1]);

%% create W - matrix of coupling (force constants)

% calculate k value for each pixel
k_2 = zeros(Nx,Ny);
G_star_mod = sqrt(Slice.^2+Slice_im.^2);
for i = 1:Nx
    for j= 1:Ny
        if (G_star_mod(i,j)==0)
            k_2(i,j) = k_0;
        else
            k_2(i,j) = (2*G_star_mod(i,j)^2)/((Slice(i,j)+G_star_mod(i,j))*ro*dx^2);
        end
    end
end

%fill W matrix (see function bellow)
W= W_matrix_create(k_2);

%% Set attenuation (see function below)

gamma_attenuation = set_attenuation(Slice,Slice_im,k_2,F0,in_border,angular_freq,ro,dx);


%% Main loop

dt = 0.001; %time step
t_end = 0.5; %simulation time

%perform simulation (see function below)
UU = CHO_simulation(angular_freq, Nx, Ny, W, gamma_attenuation, F0vec, dt, t_end);

toc %calculate elapsed time 

%% gif

% for i=1:N
%     fig = figure(200); clf
%     contourf(real(UU_8(:,:,i))); colorbar
%     title('wave propagation');
% %    caxis([min(min(real(UU_8(:,:,i))))*0.01,0.01*max(max(real(UU(:,:,i))))]);
%     drawnow
%     frame = getframe(fig);
%     im{i} = frame2im(frame);
% end

%% to save gif+file

% filename = "./Save/wave_propagation.gif"; % Specify the output file name
% for idx = 1:N
%     [L,map] = rgb2ind(im{idx},256);
%     if idx == 1
%         imwrite(L,map,filename,"gif","LoopCount",Inf,"DelayTime",2);
%     else
%         imwrite(L,map,filename,"gif","WriteMode","append","DelayTime",2);
%     end
% end

%% Fitting time series to the sine curve

u_final = fit_to_sine(UU,angular_freq,0,dt);
%u_final = fit_to_sine_fft(UU);

%%
u_final = u_final.*mask;

%% Evaluation
bias = 0; %7;

if plane == 1
    ref = double(orig_re_3+bias).*mask;
    sim = real(u_final);
else
    ref = double(orig_re_1+bias).*mask;
    sim = -real(u_final);
end

max_ref = max(max(ref));
min_ref = min(min(ref));
    
%% plot real and imaginary parts
figure(107); clf
contourf(sim,15); 
cbar = colorbar;
l = ylabel(cbar,'Displacement (μm)','Rotation',270);
l.Position(1) =3;
caxis([min(min(real(u_final))),max(max(real(u_final)))]);
xlabel('Number of pixels'); 
ylabel('Number of pixels'); 
title('Real part of displacement');

figure(109); clf
contourf(ref,15);
cbar = colorbar;
l = ylabel(cbar,'Displacement (μm)','Rotation',270);
l.Position(1) =3;
%caxis([min(min(real(u_final))),max(max(real(u_final)))]);
xlabel('Number of pixels'); 
ylabel('Number of pixels'); 
title('Reference displacement along x axis (left to right)');

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

disp(['RMSE for real part = ' num2str(RMSE_re) '; Max ref value = ' num2str(max_ref) '; Min ref value = ' num2str(min_ref)]);

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

disp(['MAE  = ' num2str(MAE_re) ', Bias = ' num2str(bias) ' MAE % = ' num2str(MAE_pers)]);

%% Plot orig vs sim plot along x (real)

figure(112); clf
hold on
scatter(1:1:Ny, squeeze(sim(50,:)));
scatter(1:1:Ny, squeeze(ref(50,:)));
%ylim([-5 10])
legend('sim','ref', 'fontsize',12)
%caxis([min(min(real(u_final))),max(max(real(u_final)))]);
xlabel('Number of pixel', 'fontsize',15); 
ylabel('Displacement (μm)','fontsize',15); 
%title('Displacement for 50th row');
hold off

%% plot simulation vs reference
sim_reshape = reshape(sim,[Nx*Ny,1]);
ref_reshape = reshape(ref,[Nx*Ny,1]);

X = [ones(length(ref_reshape),1) ref_reshape];
b2 = X\sim_reshape;
lin2 = X*b2;

figure(1010);
clf
plot (ref_reshape, sim_reshape,'o');
hold on
plot (ref_reshape, lin2);
xlabel('Reference displacement', 'fontsize',14);
ylabel('Simulated displacement', 'fontsize',14);
lgd = legend('Data','Fit');


%% SSIM
sim = mat2gray(sim);
ref = mat2gray(ref);

figure(1109);
imshow(sim);

figure(1110);
imshow(ref);

SSIM_re = ssim(sim,ref);
disp(['SSIM for real part = ' num2str(SSIM_re)]);

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
            in_border(i,j)=1; %0;
        elseif  (mask(i-1,j)==0 || mask(i+1,j)==0) && mask(i,j)==1
            in_border(i,j)=1;% 0;
        else
            in_border(i,j)=0; %mask(i,j);
        end
    end
end
%in_border = mask-in_border;
end

function F0 = set_force(scheme,in_border)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set_force creates a mask of external driving source
%
%Input: 
%   scheme - exitation scheme of the brain (1 or 2) 
%   in_border - inner boundary of the brain
%   
%Output:
%   F0 - mask 0f external force 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F0=in_border(:,:);

if scheme == 1
    %first option
    F0(:,25:end)=0;
else 
    %second option
    F0(:,25:110)=0;
    F0(:,110:end)=F0(:,110:end)*(1);
end
end

function gamma_attenuation = set_attenuation(Slice, Slice_im,k_2,F0,in_border,angular_freq, ro,dx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%gamma_attenuation creates an attenuation matrix for 2D image 
%
%Input: 
%   Slice     - matrix of real part of complex shear modulus
%   Slise_im  - matrix of imaginary part of complex shear modulus
%   k_2       - matrix with k^2 values (force constants)
%   F0        - mask of external driving source
%   in_border - mask of inner boundary of the brain
%   ro        - density of the tissues
%   dx        - distance brtween two ocsillators
%   
%Output:
%   gamma_attenuation - attenuation matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = size(k_2,1); Ny = size(k_2,2);
atten=zeros(Nx,Ny);
for i = 1:Nx
    for j=1:Ny
        G_star_mod = sqrt(Slice(i, j)^2 + Slice_im(i,j)^2);
        if Slice(i,j)==0
            %attenuation outside the brain
            atten(i,j)=2000; 
        else
            if in_border(i,j)==1 && F0(i,j)==0
                %BC on the boundaries of the brain 
                atten(i,j) = sqrt(4*k_2(i,j));
            else
                %attenuation coefficients calculated from complex shear
                %modulus
                alfa = angular_freq*sqrt(ro*1000*(G_star_mod - Slice(i,j)))/(sqrt(2*1000*G_star_mod^2));
                atten(i,j) = angular_freq*(1-exp(-alfa*dx));
            end
        end
    end
end

atten_r=reshape(atten,[Nx*Ny,1]);
%set diagonal matrix
gamma_attenuation=zeros(Nx*Ny);
for i=1:Nx*Ny
    gamma_attenuation(i,i) = atten_r(i);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fit_to_sine fits 2D time series target frequency sine curve using fft
% 
%Input: 
%   UU           : 2D time series
%   
%Output:
%   u_final      : matrix, where each pixel is a sine curve 
%                  in complex representation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fit_to_sine fits 2D time series to sime curve using fft
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

for i = 1:Nx
    for j = 1:Ny
        vect = squeeze(UU(i,j,:));
        N = length(vect);
        d_fft=fft(vect); %Fast Fourier Transform
        [~,I] = max(abs(d_fft(1:ceil(N/2)+1))); %index of max frequency 
        
        R = abs(d_fft(I)/N); %magnitude 
        phi = angle(d_fft(I)); %phase
        
        % complex representation
        u_final(i,j) = R*(cos(phi)+sqrt(-1)*sin(phi));
    end
end
end