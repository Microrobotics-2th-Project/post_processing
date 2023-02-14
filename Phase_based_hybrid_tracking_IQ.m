%% ----- IMPORT DATA -----

clear all
clc

% file ='C:\Users\Andre\Desktop\Stefano\TELEMED SW\SDK\_USGFW2 SDK 3.37\SDK\samples_cpp_vs2017\Giovanni\RF_process_app1\Experiment_file.bin'; fps=66; frot=1; fvib=8; 
file = 'Experiment_file1.bin';  fps=66; frot=1; fvib=5;


% --------DATA PREPARATION---------
[Data,parameters,tranducer,success]=import_RFdata2MATLAB_pane_mod(file);
% [Data]=import_RFdata2MATLAB_pane_mod2(file);
Fs=40000000;
Ftx=10000000;
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
    Au=Data{1,j};
    I(:,:,j)=Au(1,:,:);
    Q(:,:,j)=Au(2,:,:);
end

% I=[];
% % RF data
% for j=1:length(Data)
%     I(:,:,j)=Data{1,j};
% end

% % --ROI selection--
% I=I(100:1000,30:160,:);
% 
% ph2 = I;

% p=hilbert(I);              % analytic signal (along y)
% ph=unwrap(angle(p),[],3);  % acoustic phase signal (istantaneous phase)
% fd=diff(ph,1,3)*fps;       % differential phase signal (istantanous frequency)

ph2=unwrap(atan2(Q,I),[],3);
% ph2=atan2(Q,I);
fd2=diff(ph2,1,3)*fps; 

T=round(fps);   % magnetic signal start
T0=T;

PRF=fps;
max_w=PRF/2*2*pi*const; % m/s


%----MICROROBOT PARAMETERS----
r_eff=5.5e-4; % expected MR radius
tp_size=6*r_eff;  % template size
%---TEMPLATE DEFINITION---
square=0;
circ=not(square);
tp=zeros(round(tp_size/res_y),round(tp_size/res_x));
x=1:1:size(tp,2);
y=1:1:size(tp,1);
[X,Y]=meshgrid(x,y);
vy=sqrt((X-size(tp,2)/2).^2+(Y-size(tp,1)/2).^2).*cos(atan2(Y-size(tp,1)/2,X-size(tp,2)/2));
%vy = -(2*(X-size(tp,2)/2).^2 + (1/10)*(Y-size(tp,1)/2).^2)+100;
if(circ)
    vy((((X-size(tp,2)/2).^2)+(((Y-size(tp,1)/2)*res_y/res_x).^2))>(round(r_eff/res_x))^2)=0;
end
if(square)
    dx=round(2*r_eff/res_x);
    dy=round(2*r_eff/res_y);
    vy(1:(y(end)-dy)/2,:)=0;
    vy(:,1:(x(end)-dx)/2)=0;
    vy(:,(x(end)+dx)/2:end)=0;
    vy((y(end)+dy)/2:end,:)=0;
end
tp=vy;

% %------PLOT TEMPLATE------
% figure
% imagesc(vy)
% xv=5:5:20;
% yv=10:10:80;
% xticklabels(150e-3*(xv));
% xticks(xv)
% yticklabels(37.5e-3*((yv)));
% yticks(yv)
% %-------------------------
% I = ph2;