%% ----- IMPORT DATA -----

clear all
clc

%------------------OLD DATA--------------------
% -----------9GIU---------------
% file='14.18.37_09-06-2020_L15-7H40-A5.bin'; fps=119;fm=3; % no field
% file='15.38.41_08-06-2020_L15-7H40-A5.bin'; fps=119;fm=3; % vibrating MR
% file='17.27.53_02-12-2019_L15-7H40-A5.bin'; fps=138.2; fm=3;  % MR 2
% file='18.43.35_09-11-2020_L15-7H40-A5.bin'; fps=102; fm=2; % vibrating in place 
% file='13.36.59_10-06-2020_L15-7H40-A5.bin'; fps=119;fm=3; % tracking 0.1
% file='13.42.14_10-06-2020_L15-7H40-A5.bin'; fps=119;fm=3; % tracking 0.5
% file='13.43.44_10-06-2020_L15-7H40-A5.bin'; fps=119;fm=3; % tracking 1

%------------------ROTATION IN PLACE--------------------
% -----------22OTT---------------
% ------3mm-----
% file='14.49.43_22-10-2020_L15-7H40-A5.bin'; fps=104;frot=0.3; % ROT 3mm 0.3Hz (0.18)
% file='14.52.03_22-10-2020_L15-7H40-A5.bin'; fps=104;frot=0.6; fvib = 1;% ROT 3mm 0.6Hz (0.26)
% file='14.55.08_22-10-2020_L15-7H40-A5.bin'; fps=104;frot=0.9; % ROT 3mm 0.9Hz (0.37)
% file='15.00.52_22-10-2020_L15-7H40-A5.bin'; fps=104;fm=1.2;r_eff=1.5e-3; % ROT 3mm 1.2Hz (1.6e-3 m/s)
% file='15.02.06_22-10-2020_L15-7H40-A5.bin'; fps=104;fm=1.5;r_eff=1.5e-3; % ROT 3mm 1.5Hz
% ------1mm-----
% file='15.04.30_22-10-2020_L15-7H40-A5.bin'; fps=104;frot=0.3; % ROT 1mm 0.3Hz (0.15)
% file='15.05.46_22-10-2020_L15-7H40-A5.bin'; fps=104;frot=0.6; % ROT 1mm 0.6Hz (0.3)
% file='15.06.56_22-10-2020_L15-7H40-A5.bin'; fps=104;frot=0.9; % ROT 1mm 0.9Hz (0.35)
% file='15.08.28_22-10-2020_L15-7H40-A5.bin'; fps=104;fm=1.2; % ROT 1mm 1.2Hz (8e-4 m/s)
% file='15.09.36_22-10-2020_L15-7H40-A5.bin'; fps=104;fm=1.5; % ROT 1mm 1.5Hz (8e-4 m/s)
% -----------29OTT---------------
% ------3mm ROT-----
% file='18.39.08_29-10-2020_L15-7H40-A5.bin'; fps=104;frot=0.2;fvib=1; % 0.2Hz (0.14)
% file='18.42.53_29-10-2020_L15-7H40-A5.bin'; fps=104;frot=0.4;fvib=1; % 0.4Hz (0.28)
% file='18.45.03_29-10-2020_L15-7H40-A5.bin'; fps=104;frot=0.6;fvib=1; % 0.6Hz (0.4)
% file='18.46.53_29-10-2020_L15-7H40-A5.bin'; fps=104;fm=0.8; % 0.8Hz (B)
% file='18.48.38_29-10-2020_L15-7H40-A5.bin'; fps=104;fm=1; % 1Hz (B)
% ------1mm ROT-----
% file='18.41.07_29-10-2020_L15-7H40-A5.bin'; fps=104;frot=0.2;fvib=1; % 0.2Hz (0.14) 
% file='18.43.27_29-10-2020_L15-7H40-A5.bin'; fps=104;frot=0.4;fvib=1; % 0.4Hz (0.28)
% file='18.45.33_29-10-2020_L15-7H40-A5.bin'; fps=104;frot=0.6;fvib=1; % 0.6Hz (0.4)
% file='18.47.26_29-10-2020_L15-7H40-A5.bin'; fps=104;fm=0.8; % 0.8Hz (B)
% file='18.49.13_29-10-2020_L15-7H40-A5.bin'; fps=104;fm=1; % 1Hz (B)


%------------------ROLLING LOCOMOTION------------------------
%-----15cm------
% file='17.37.21_28-09-2020_L15-7H40-A5.bin'; fps=103;fm=0.5; % roll 0.5 Hz
% file='17.40.28_28-09-2020_L15-7H40-A5.bin'; fps=103;fm=0.5; % roll 0.5 Hz
% file='17.42.49_28-09-2020_L15-7H40-A5.bin'; fps=103;fm=1; % roll 1 Hz
% file='17.43.09_28-09-2020_L15-7H40-A5.bin'; fps=103;fm=1; % roll 1 Hz
% file='17.44.53_28-09-2020_L15-7H40-A5.bin'; fps=103;fm=1.5; % roll 1.5 Hz
% file='17.45.30_28-09-2020_L15-7H40-A5.bin'; fps=103;fm=1.5; % roll 1.5 Hz
%-----20cm------
% file='17.48.11_28-09-2020_L15-7H40-A5.bin'; fps=103;fm=0.5; % in vessel 0.5 Hz (5e-4 m/s)
% file='17.49.11_28-09-2020_L15-7H40-A5.bin'; fps=103;fm=0.5; % in vessel 0.5 Hz (5e-4 m/s)
% file='17.50.18_28-09-2020_L15-7H40-A5.bin'; fps=103;fm=1; % in vessel 1 Hz (7e-4 m/s)
% file='17.50.35_28-09-2020_L15-7H40-A5.bin'; fps=103;fm=1; % in vessel 1 Hz (7e-4)
% file='17.52.27_28-09-2020_L15-7H40-A5.bin'; fps=103;fm=1.5; % in vessel 1.5 Hz
% file='17.52.41_28-09-2020_L15-7H40-A5.bin'; fps=103;fm=1.5; % in vessel 1.5 Hz
% -----------30SET---------------
% file='17.31.13_30-09-2020_L15-7H40-A5.bin'; fps=103;fm=1; % 10cm low
% file='17.32.21_30-09-2020_L15-7H40-A5.bin'; fps=103;fm=1; % 10cm high
% file='16.33.29_30-09-2020_L15-7H40-A5.bin'; fps=103;fm=1; % free flow
% file='17.25.50_30-09-2020_L15-7H40-A5.bin'; fps=103;fm=1; % 5cm 1.5V
% file='17.26.32_30-09-2020_L15-7H40-A5.bin'; fps=103;fm=1; % 5cm 2V
% file='17.27.24_30-09-2020_L15-7H40-A5.bin'; fps=103;fm=1; % 5cm 2.5V
% file='17.28.03_30-09-2020_L15-7H40-A5.bin'; fps=103;fm=1; % 5cm 3V
% -----------31SET---------------
% file='13.03.12_01-10-2020_L15-7H40-A5.bin'; fps=103;fm=1; % 10cm 1V
% file='13.04.16_01-10-2020_L15-7H40-A5.bin'; fps=103;fm=1; % 10cm 1.5V
% file='13.05.37_01-10-2020_L15-7H40-A5.bin'; fps=103;fm=1; % 10cm 2V
% file='13.08.38_01-10-2020_L15-7H40-A5.bin'; fps=103;fm=1;r_eff=5e-4; % 10cm 2.5V Good
% -----------16OTT---------------
% file='16.36.09_16-10-2020_L15-7H40-A5.bin'; fps=103;fm=1; % MR no stream (Good)
% file='16.36.25_16-10-2020_L15-7H40-A5.bin'; fps=103;fm=1; % MR no stream2
% file='16.16.13_16-10-2020_L15-7H40-A5.bin'; fps=103;fm=1; % MR strong stream
% file='16.37.20_16-10-2020_L15-7H40-A5.bin'; fps=103;fm=1; % MR strong stream2
% file='16.38.14_16-10-2020_L15-7H40-A5.bin'; fps=103;fm=1; % MR 1.8V stream
% file='16.38.45_16-10-2020_L15-7H40-A5.bin'; fps=103;fm=1; % MR 1.8V stream2
% file='16.39.31_16-10-2020_L15-7H40-A5.bin'; fps=103;fm=1; % MR 1.8V stream3 (TOP)
% -----------22OTT---------------
% ------coated glycerol-----
% file='16.41.13_22-10-2020_L15-7H40-A5.bin'; fps=103;frot=1.5;fvib=1; % low flux      Bmode1(high contrast)
% file='16.41.48_22-10-2020_L15-7H40-A5.bin'; fps=103;fm=1.5;fvib=1; % high flux
% file='16.42.20_22-10-2020_L15-7H40-A5.bin'; fps=103;fm=1.5; % mid flux
% --------non-coated water------
% file='16.29.52_22-10-2020_L15-7H40-A5.bin'; fps=103;fm=1.5;fvib=1;
% ------non-coated glycerol-----
file='16.49.29_22-10-2020_L15-7H40-A5.bin'; fps=103;frot=1.5;fvib=1; % no-flux (G)     Bmode5
% file='16.52.39_22-10-2020_L15-7H40-A5.bin'; fps=103;frot=1.5;fvib=1; % flux (G) 
% file='16.53.33_22-10-2020_L15-7H40-A5.bin'; fps=103;frot=1.5;fvib=1; % mid flux (VG) --
% file='16.54.11_22-10-2020_L15-7H40-A5.bin'; fps=103;frot=1.5;fvib=1; % high flux (VG)
% file='16.55.08_22-10-2020_L15-7H40-A5.bin'; fps=103;frot=1.5;fvib=1; % high flux (G)
% -----------29OTT---------------
% ------coated glycerol-----
% file='14.44.51_29-10-2020_L15-7H40-A5.bin'; fps=103;frot=1;fvib=1; % no flux 
% file='14.49.27_29-10-2020_L15-7H40-A5.bin'; fps=103;frot=1;fvib=1; % 1.2V flux          (good)
% file='14.46.18_29-10-2020_L15-7H40-A5.bin'; fps=103;frot=1;fvib=1; % 2V flux (B)
% ------non-coated glycerol-----
% file='14.55.17_29-10-2020_L15-7H40-A5.bin'; fps=103;frot=1; fvib=1;% noflux (G) --    (Bmode4)
% file='14.57.05_29-10-2020_L15-7H40-A5.bin'; fps=103;frot=1; fvib=1; % 1.2V flux 

%------------------ROLLING + VIBRATIONS------------------------
% -----------19NOV---------------
% file='16.54.28_19-11-2020_L15-7H40-A5.bin'; fps=119; frot=1; fvib=8.3; % wid 3 water
% file='17.05.13_19-11-2020_L15-7H40-A5.bin'; fps=119; frot=1; fvib=8.3; % wid 3 glyc
% file='16.36.41_19-11-2020_L15-7H40-A5.bin'; fps=119; frot=1; fvib=4.95; % wid 5 water
% file='17.12.06_19-11-2020_L15-7H40-A5.bin'; fps=119; frot=1; fvib=4.95; % wid 5 glyc
% file='16.44.17_19-11-2020_L15-7H40-A5.bin'; fps=119; frot=1; fvib=2.45; % wid 10 water
% file='17.16.03_19-11-2020_L15-7H40-A5.bin'; fps=119; frot=1; fvib=2.45; % wid 10 glyc
% file='17.24.10_19-11-2020_L15-7H40-A5.bin'; fps=119; frot=1; fvib=4.95; % combined wid5 glyc (B)
% file='17.29.59_19-11-2020_L15-7H40-A5.bin'; fps=119; frot=1; fvib=4.95; % combined wid5 glyc
% file='17.39.47_19-11-2020_L15-7H40-A5.bin'; fps=119; frot=1; fvib=4.95; % combined wid5 glyc (G) --  Bmode3(possible)
% file='17.22.15_19-11-2020_L15-7H40-A5.bin'; fps=119; frot=1; fvib=4.95; % combined no focus
% -----------23NOV---------------
% file='17.14.06_23-11-2020_L15-7H40-A5.bin'; fps=119; frot=1; fvib=4.95; % wid 5 flow
% file='17.17.01_23-11-2020_L15-7H40-A5.bin'; fps=119; frot=1; fvib=4.95; % wid 5 flow (G) --
% file='17.19.01_23-11-2020_L15-7H40-A5.bin'; fps=119; frot=1; fvib=8.3; % wid 3 flow (mmh)
% file='17.19.53_23-11-2020_L15-7H40-A5.bin'; fps=119; frot=1; fvib=8.3; % wid 3 flow (mmh)


% --------DATA PREPARATION---------
tic;
[Data,parameters,tranducer,success]=import_RFdata2MATLAB_pane(file);
toc
%Data = import_RFdata2MATLAB_Test('Test_Results.bin');
Fs=1e9/parameters{1,250}.Sampling_period_ns;
Ftx=parameters{1,250}.tx_frequency;
% fps=parameters{1,250}.frame_rate;
% fps=119;
c=1485;
% Ftx=7500000;
% c=1500;
lambda=c/(Ftx);
const=lambda/(4*pi);
const=const*1.4;
res_x=150e-6;
res_y=37e-6;
% const1=3.7e-9/0.8e-3;


I=[];
% RF data
for j=1:length(Data)
    I(:,:,j)=Data{1,j};
end

% % --ROI selection--
% I=I(100:1000,30:160,:);
% [M,N,Q] = size(I);

p=hilbert(I);              % analytic signal (along y)
%write p in binary file
[M,N,Q] = size(p);
ens = 34;
index = 1:43;
% to_write = [fps,fvib,M,N];
% %for i=1:43  %iteration
%     i = 21
%     shift = 1+ens*i;
%     for j=0:ens             %cineloop
%         for m=1:M           %row
%             for n=1:N       %column
%                 re = real(p(m,n,shift+j));
%                 im = imag(p(m,n,shift+j));
%                 to_write = [to_write,re,im];
%             end
%         end
%     end
% %end
ph=unwrap(angle(p),[],3);  % acoustic phase signal (istantaneous phase)
fd=diff(ph,1,3)*fps;       % differential phase signal (istantanous frequency)
[M,N,Q] = size(ph);
% to_write2 = [fps,fvib,M,N];
% to_write = [M,N];

T=round(fps);   % magnetic signal start
T0=T;

PRF=fps;
max_w=PRF/2*2*pi*const; % m/s