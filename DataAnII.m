%Test to wrap together different cineloop images
clear all
close all
myfig

%% Load RF Data
file = 'Experiment_file4_flux_small.bin';
% file = 'Experiment_file1.bin';
[Data,parameters,tranducer,success]=import_RFdata2MATLAB_pane_mod(file);
Fs=40000000;
Ftx=10000000;
c=1485;
lambda=c/(Ftx);
const=lambda/(4*pi);
const=const*1.4;
res_x=150e-6;
res_y=37e-6;
I=[];
% RF data
for j=1:length(Data)
    Au=Data{1,j};
    I(:,:,j)=Au(1,:,:);
    Q(:,:,j)=Au(2,:,:);
end

%% Load Robot Pose data
%pixel to millimeter relationship:
pixel_height_RF = 0.038;                        %[mm]
pixel_width = 0.150;                            %[mm]
pixel_height_B = 4*pixel_height_RF;             %[mm]
load('exp4_flux_small.mat');
% load('exp1.mat');
indexes = [];
for m = 2:length(xMR)
    if((xMR(m)~=xMR(m-1))||(zMR(m)~=zMR(m-1)))
        indexes = [indexes,m];
    end
end
% indexes = indexes(4:end);
x_img = pose(indexes,2)+50;     %for exp2 30/06 use +90, for exp1 use +50                   
y_img = pose(indexes,3)-210;
x_orig = [];
y_orig = [];
for k=length(x_img):-1:2
    x_orig = [x_orig, round(-(x_img(k)-x_img(k-1))/pixel_width)];
    y_orig = [y_orig, round((y_img(k)-y_img(k-1))/(pixel_height_RF))];
end
%probably first 10-15 cineloop are to throw away
% x_img = x_img(3:end-3);
% y_img = y_img(3:end-3);
% %for exp1 30/06
% x_img = x_img(30:end);
% y_img = y_img(30:end);

%% Compress Cineloops
[M,N,P] = size(I);
cineloop_dimension = 60;
cineloop_nr = P/cineloop_dimension;
% cineloop_nr = cineloop_nr -29;  %for exp1 30/06
if(cineloop_nr>length(indexes))
    cineloop_nr = length(indexes);
end
meanI = zeros(M,N,cineloop_nr);
%meanI = [];
averageI = zeros(M,N);
%averageI = [];
for j = 1:cineloop_nr
    for i = 1+(j-1)*cineloop_dimension:j*cineloop_dimension
%         if(i==1+(j-1)*cineloop_dimension)
%             averageI = rf2bmode(I(:,:,i));
%         else
%             averageI = averageI + rf2bmode(abs(I(:,:,i)+Q(:,:,i)));
%         end
    averageI = averageI + I(:,:,i);
    end
    meanI(:,:,j) = averageI./cineloop_dimension;
    averageI = zeros(M,N);
end

%% Movement Image Reconstruction
%Method 1: merge rf2bmode images of cineloop mean frame
A = rf2bmode(meanI(:,:,1));
Reconstruction = zeros(size(A));
start_x = 1;
start_y = 1;
start_y2 = 1;
%Exp 1 11_05 comment if(k<10||k>cineloop_nr-40) part
%for exp1 30/06 start 24 end 0
%for exp2 30/06 start 0 end 0
%for exp3 30/06 start 105 end 0
%for exp4 30/06 start 9 end 0
%for exp5 30/06 start 26 end 0
%for exp6 30/06 start 100 end 0
offset_start = 9;
offset_end = 0;
th = 15;
vessel_x = [];
vessel_y = [];
ref_y_vessel = [];
ref_x_vessel = [];
for k = 1+offset_end:cineloop_nr-offset_start
    A = rf2bmode(meanI(:,:,cineloop_nr-k+1));
    [M,N] = size(A);
    A(1:3,:) = 0;
    A(M-2:M,:) = 0;
%     A = NSRFilters(A,'ideal',40);           %Fourier Ideal Filter with 40% cutoff frequency
%     %vessel segmentation in single frame
%         A_filtered = sobel(A,th);
%         A_filtered(1:10,:) = 0;
%         A_filtered(M-5:M,:) = 0;
%         A_filtered(:,1:3) = 0;
%         A_filtered(:,N-3:N) = 0;
%         A_filtered = double(A_filtered>0);
%         for i=1:N
%             [pks,locs] = findpeaks(A_filtered(:,i));
%             if(length(locs)>1)
%                 A_filtered(locs(2):M,i) = 0;
%             end
%         end
        filter = [zeros(1,60); ones(1,60); zeros(1,60)];
        A_fil = imfilter(A,filter);
        for i = 1:N
            [pks,locs] = findpeaks(A_fil(:,i));
            vessel_x = [vessel_x;i];
            y = locs(pks==max(pks));
            pks(pks==max(pks)) = 0;
%             vessel_y = [vessel_y;(y+locs(pks==max(pks)))/2];
            vessel_y = [vessel_y;y];
        end
        ref_y_vessel = [ref_y_vessel,N];
        ref_x_vessel = [ref_x_vessel,N];
%         [vessel_yA,vessel_xA] = find(A_filtered==1);
%         ref_y_vessel = [ref_y_vessel,length(vessel_yA)];
%         ref_x_vessel = [ref_x_vessel,length(vessel_xA)];
%         vessel_x = [vessel_x;vessel_xA];
%         vessel_y = [vessel_y;vessel_yA];
    %Merge of frames
    if(k>1)
        start_x = start_x + x_orig(k-1);
        start_y = start_y + y_orig(k-1);
        start_y2 = start_y2 + floor(y_orig(k-1)/4);
        if(start_y2==0)
            start_y2=1;
        end
    end
%     if(k<10||k>cineloop_nr-40)
%         A = zeros(M,N);
%     end
    if((start_x<=0)||(start_y2<=0))
        [Mr,Nr] = size(Reconstruction);
        Ma = Mr;
        Na = Nr;
        start_x_Re = 1;
        start_y_Re = 1;
        if(start_x<=0)
            Na = Nr-start_x+1;
            start_x_Re = abs(start_x)+1;
            start_x = 1;
        end
        if(start_y2<=0)
            Ma = Mr-start_y2+1;
            start_y_Re = abs(start_y2)+1;
            start_y2 = 1;
        end
        Aux = Reconstruction;
        Reconstruction = zeros(Ma,Na);
        Reconstruction(start_y2:start_y2+M-1,start_x:start_x+N-1) = A;
        Reconstruction(start_y_Re:start_y_Re+Mr-1,start_x_Re:start_x_Re+Nr-1) = Aux;
        
%         [Ma,Na] = size(Reconstruction);
%         Aux = zeros(Ma+abs(start_x)+1,Na+abs(start_y2)+1);
%         Aux(abs(start_y2)+1:abs(start_y2)+Ma,abs(start_x)+1:abs(start_x)+Na) = Reconstruction;
%         Reconstruction = Aux;
%         Reconstruction(1:M,1:N) = A;
    else
        Reconstruction(start_y2:start_y2+M-1,start_x:start_x+N-1) = A;
    end
end
%Method 2: rf2bmode of merge of cineloop mean frame
[M,N,P] = size(I);
Reconstruction2 = zeros(M,N);
start_x = 1;
start_y = 1;
ref_coord_row = zeros(1,cineloop_nr);
ref_coord_column = zeros(1,cineloop_nr);
for k = 1+offset_end:cineloop_nr-offset_start
    if(k>1)
        start_x = start_x + x_orig(k-1);
        start_y = start_y + y_orig(k-1);
    end
    ref_coord_column(cineloop_nr-k+1) = start_x;
    ref_coord_row(cineloop_nr-k+1) = start_y;
    if((start_x<=0)||(start_y<=0))
        [Mr,Nr] = size(Reconstruction2);
        Ma = Mr;
        Na = Nr;
        start_x_Re = 1;
        start_y_Re = 1;
        if(start_x<=0)
            Na = Nr-start_x+1;
            start_x_Re = abs(start_x)+1;
            start_x = 1;
            ref_coord_column(cineloop_nr-k+1) = start_x;
            ref_coord_column(cineloop_nr-k+2:end) = ref_coord_column(cineloop_nr-k+2:end)+start_x_Re;
        end
        if(start_y<=0)
            Ma = Mr-start_y+1;
            start_y_Re = abs(start_y)+1;
            start_y = 1;
            ref_coord_row(cineloop_nr-k+1) = start_y;
            ref_coord_row(cineloop_nr-k+2:end) = ref_coord_row(cineloop_nr-k+2:end)+start_y_Re;
        end
        Aux = Reconstruction2;
        Reconstruction2 = zeros(Ma,Na);
        Reconstruction2(start_y:start_y+M-1,start_x:start_x+N-1) = meanI(:,:,cineloop_nr-k+1);
        Reconstruction2(start_y_Re:start_y_Re+Mr-1,start_x_Re:start_x_Re+Nr-1) = Aux;
%         [Ma,Na] = size(Reconstruction2);
%         Aux = zeros(Ma+abs(start_x)+1,Na+abs(start_y)+1);
%         Aux(abs(start_y)+1:abs(start_y)+Ma,abs(start_x)+1:abs(start_x)+Na) = Reconstruction2;
%         Reconstruction2 = Aux;
%         Reconstruction2(1:M,1:N) = meanI(:,:,cineloop_nr-k+1);
%         if(k<10||k>cineloop_nr-40)
%             Reconstruction2(start_y:start_y+M-1,start_x:start_x+N-1) = zeros(M,N); 
%         end
    else
        Reconstruction2(start_y:start_y+M-1,start_x:start_x+N-1) = meanI(:,:,cineloop_nr-k+1);
%         if(k<10||k>cineloop_nr-40)
%             Reconstruction2(start_y:start_y+M-1,start_x:start_x+N-1) = zeros(M,N); 
%         end
    end
end
Reconstruction2 = rf2bmode(Reconstruction2);

%% Trajectory reconstruction
%rf2bmode compress the 4 rows into 1, so we need to change y representation
%for microrobot
%first retrieve MR coordinates in cineloop frames reference
% %for exp3 30/06 and exp1 11/05
% xMR_new = round(xMR(indexes).*(1000/pixel_width) + N/2);
% yMR_new = round(zMR(indexes).*(1000/pixel_height_RF) + M/2);
% % for exp1 30/06
% xMR_new = round(xMR(indexes(30:end)).*(1000/pixel_width) + N/2);
% yMR_new = round(zMR(indexes(30:end)).*(1000/pixel_height_RF) + M/2);
% %for exp2 30/06
% xMR_new = round(xMR(indexes(3:end-3)).*(1000/pixel_width) + N/2);
% yMR_new = round(zMR(indexes(3:end-3)).*(1000/pixel_height_RF) + M/2);
%for exp4 30/06
xMR_new = round(xMR(indexes(10:end-0)).*(1000/pixel_width) + N/2);
yMR_new = round(zMR(indexes(10:end-0)).*(1000/pixel_height_RF) + M/2);
% %for exp5 30/06, segment and use only the first part
% xMR_new = round(xMR(indexes(3:end-5)).*(1000/pixel_width) + N/2);
% yMR_new = round(zMR(indexes(3:end-5)).*(1000/pixel_height_RF) + M/2);
% %for exp6 30/06
% xMR_new = round(xMR(indexes(3:end-0)).*(1000/pixel_width) + N/2);
% yMR_new = round(zMR(indexes(3:end-0)).*(1000/pixel_height_RF) + M/2);
%compress yMR
yMR_compressed = zeros(1,length(yMR_new));
yMR_plot2 = yMR_new+ref_coord_row';
for i=1:length(yMR_new)
    yMR_compressed(i) = floor(yMR_new(i)/4);
    ref_coord_row(i) = floor(ref_coord_row(i)/4);
end
%using ref_coord we can obtain the correct coordinates for MR in
%Reconstruction
xMR_plot = xMR_new + ref_coord_column';
yMR_plot = yMR_compressed+ref_coord_row;

% %fit for exp1 11/05
% xMR_plot_ref = mean(xMR_new(14:end-1))+ref_coord_column';
% 
% error_Vis_x = abs((N/2)-mean(xMR_new(14:end-1))).*pixel_width;
% error_Vis_y = abs((floor(M/4)/2)-mean(yMR_compressed(14:end-1))).*pixel_height_B;
% 
% %results from fit of segmented vessel
% offset = 35.32;
% coeff1 = 0.1887;
% coeff2 = -0.0002457;
% coeff3 = 1.493e-07;
% yMR_plot_ref2 = coeff3*(xMR_plot_ref.^3) +coeff2*(xMR_plot_ref.^2) + coeff1*xMR_plot_ref + offset;    %segmented vessel
% yMR_plot_ref = (yMR_plot_ref2-6)';     %align with MR

% %fit for exp1 30/06
% xMR_plot_ref = mean(xMR_new(25:end))+ref_coord_column';
% 
% error_Vis_x = abs((N/2)-mean(xMR_new(25:end))).*pixel_width;
% error_Vis_y = abs((floor(M/4)/2)-mean(yMR_compressed(25:end))).*pixel_height_B;
% 
% %results from fit of segmented vessel
% offset = 46+7;
% coeff = -0.02159;
% yMR_plot_ref2 = coeff*xMR_plot_ref + offset;    %segmented vessel
% yMR_plot_ref = (yMR_plot_ref2-7)';     %align with MR

% %fit for exp2 30/06
% xMR_plot_ref = mean(xMR_new(20:end-10))+ref_coord_column';
% 
% %results from fit of segmented vessel
% offset = 71;
% coeff = -0.03772;
% yMR_plot_ref2 = coeff*xMR_plot_ref + offset;    %segmented vessel
% yMR_plot_ref = (yMR_plot_ref2-7)';     %align with MR

% %fit for exp3 30/06
% xMR_plot_ref = mean(xMR_new(106:end))+ref_coord_column';
% 
% %results from fit of segmented vessel
% offset = 48;
% coeff = -0.02394;
% yMR_plot_ref2 = coeff*xMR_plot_ref + offset;    %segmented vessel
% yMR_plot_ref = (yMR_plot_ref2-8)';     %align with MR

%fit for exp4 30/06
values = polyfit(ref_coord_column(10:end),ref_coord_row(10:end),1);
ref_fit = values(1)*ref_coord_column(10:end)+values(2);
ref_fit = [zeros(1,9),ref_fit];

xMR_plot_ref = mean(xMR_new(10:end))+ref_coord_column';
yMR_plot_ref = mean(yMR_compressed(10:end))+ref_fit;

%results from fit of segmented vessel
offset = 56;
coeff = -0.0393;
yMR_plot_ref2 = coeff*xMR_plot_ref + offset;    %segmented vessel
yMR_plot_ref = (yMR_plot_ref2-7)';     %align with MR

% %fit for exp5 30/06
% values = polyfit(ref_coord_column(28:end),ref_coord_row(28:end),1);
% ref_fit = values(1)*ref_coord_column(28:end)+values(2);
% ref_fit = [zeros(1,27),ref_fit];
% 
% xMR_plot_ref = mean(xMR_new(28:end))+ref_coord_column';
% yMR_plot_ref = mean(yMR_compressed(28:end))+ref_fit;
% 
% %results from fit of segmented vessel
% offset = 40;
% coeff = -0.01194;
% yMR_plot_ref2 = coeff*xMR_plot_ref + offset;    %segmented vessel
% yMR_plot_ref = (yMR_plot_ref2-4)';     %align with MR

% %fit for exp6 30/06
% values = polyfit(ref_coord_column(101:end-5),ref_coord_row(101:end-5),1);
% values(1) = values(1) - 0.01;
% ref_fit = values(1)*ref_coord_column(101:end-5)+values(2);
% ref_fit = [zeros(1,100),ref_fit,zeros(1,5)];
% 
% xMR_plot_ref = mean(xMR_new(101:end-5))+ref_coord_column';
% yMR_plot_ref = mean(yMR_compressed(101:end-5))+ref_fit;

% %results from fit of segmented vessel
% offset = 58.49;
% coeff = -0.03443;
% yMR_plot_ref2 = coeff*xMR_plot_ref + offset;    %segmented vessel
% yMR_plot_ref = (yMR_plot_ref2-6)';     %align with MR

%% Load Offline results 11/05
load('Exp1_Offline.mat')
offline_xMR = c_x(8:end-4)+ref_coord_column(14:end);
offline_yMR = c_y(8:end-4);
for i=1:length(offline_yMR)
    offline_yMR(i) = floor(offline_yMR(i)/4);
end
offline_yMR = offline_yMR+ref_coord_row(14:end);

%% Plots
figure
subplot(2,1,1)
plot(x_img,y_img)
grid on
xlim([0,140])
ylim([0,20])
ylabel('[mm]')
xlabel('[mm]')
daspect([1 1 1])
title('MR Trajectory reconstruction from robot pose')
subplot(2,1,2)
m = mean(yMR_plot2);
for i=1:length(yMR_plot2)
    yMR_plot2(i) = m+sign(m-yMR_plot2(i))*abs(m-yMR_plot2(i));
end
plot(xMR_plot(14:end-1).*pixel_width,yMR_plot2(14:end-1).*(pixel_height_RF))    %exp1 15/05
% plot(xMR_plot(25:end-4).*pixel_width,yMR_plot2(25:end-4).*(pixel_height_RF))    %exp2 30/06
% plot(xMR_plot.*pixel_width,yMR_plot2.*(pixel_height_RF))                        %exp1 30/06 exp3 exp5
% plot(xMR_plot(10:end).*pixel_width,yMR_plot2(10:end).*(pixel_height_RF))          %exp4 30/06
%for exp3 30/06 segment data in the different parts and use only the first
%locomotion section
grid on
xlim([0,140])
ylim([0,20])
ylabel('[mm]')
xlabel('[mm]')
daspect([1 1 1])
title('MR Trajectory reconstruction from xMR, zMR')

figure
plot(1+(100-x_img),y_img+2,'b')
hold on
plot(1+(100-xMR_plot(14:end-1).*pixel_width),yMR_plot2(14:end-1).*(pixel_height_RF),'r')    %exp1 15/05
grid on
xlim([0,120])
ylim([0,20])
ylabel('y [mm]')
xlabel('x [mm]')
legend('ee pose','MR trajectory')
daspect([1 1 1])


figure
imagesc(Reconstruction);
colormap('gray');
hold on
% plot(xMR_plot(14:end-7),yMR_plot(14:end-7),'y') %Exp1 11_05
% plot(xMR_plot_ref(14:end-7),yMR_plot_ref(14:end-7),'g')         %Exp1 11/05
% plot(xMR_plot_ref(14:end-7),yMR_plot_ref2(14:end-7),'r')         %Exp1 11/05
% plot(offline_xMR, offline_yMR,'c')
% plot(xMR_plot(25:end-10),yMR_plot(25:end-10),'y') %Exp2 30/06
% plot(xMR_plot_ref(25:end-10),yMR_plot_ref(25:end-10),'g')         %Exp2 30/06
% plot(xMR_plot_ref(25:end-10),yMR_plot_ref2(25:end-10),'r')         %Exp2 30/06
% plot(xMR_plot(106:end),yMR_plot(106:end),'y')                         %Exp3 30/06
% plot(xMR_plot_ref(106:end),yMR_plot_ref(106:end),'g')         %Exp3 30/06
% plot(xMR_plot_ref(106:end),yMR_plot_ref2(106:end),'r')         %Exp3 30/06
% plot(xMR_plot(25:end),yMR_plot(25:end),'y')                         %Exp1 30/06
% plot(xMR_plot_ref(25:end),yMR_plot_ref(25:end),'g')         %Exp1 30/06
% plot(xMR_plot_ref(25:end),yMR_plot_ref2(25:end),'r')         %Exp1 30/06
% plot(xMR_plot(101:end-5),yMR_plot(101:end-5),'y')   %Exp 6 30/06
% plot(xMR_plot_ref(101:end-5),yMR_plot_ref(101:end-5),'g')         %Exp6 30/06
% plot(xMR_plot_ref(101:end-5),yMR_plot_ref2(101:end-5),'r')         %Exp6 30/06
% plot(xMR_plot(28:end),yMR_plot(28:end),'y')           %exp5
% plot(xMR_plot_ref(28:end),yMR_plot_ref(28:end),'g')         %Exp5 30/06
% plot(xMR_plot_ref(28:end),yMR_plot_ref2(28:end),'r')         %Exp5 30/06
plot(xMR_plot(10:end),yMR_plot(10:end),'y')         %Exp4 30/06
plot(xMR_plot_ref(10:end),yMR_plot_ref(10:end),'g')         %Exp4 30/06
plot(xMR_plot_ref(10:end),yMR_plot_ref2(10:end),'r')         %Exp4 30/06
[Mr,Nr] = size(Reconstruction);
height = round(Mr*pixel_height_B);
width = round(Nr*pixel_width);
pbaspect([width height 1])
title('Reconstruction of travelled path. Method 1')
xlabel([num2str(width),' [mm]'])
ylabel([num2str(height),' [mm]'])
legend('Tracked MR','Ground truth','Segmented vessel','Offline')
xticklabels(num2cell(round(str2double(xticklabels).*pixel_width)))
yticklabels(num2cell(round(str2double(yticklabels).*pixel_height_B)))


figure
imagesc(Reconstruction2);
colormap('gray');
hold on
% plot(xMR_plot(14:end-1),yMR_plot(14:end-1),'y') %Exp1 11_05
% plot(xMR_plot(21:end-4),yMR_plot(21:end-4),'y') %Exp2 30/06
% plot(xMR_plot,yMR_plot,'y')                         %Exp1 30/06 exp3
% plot(xMR_plot(27:end),yMR_plot(27:end),'y')           %exp5
plot(xMR_plot(10:end),yMR_plot(10:end),'y')         %Exp4 30/06
[Mr,Nr] = size(Reconstruction2);
height = round(Mr*pixel_height_B);
width = round(Nr*pixel_width);
pbaspect([width height 1])
title('Reconstruction of travelled path. Method 2')
xlabel([num2str(width),' [mm]'])
ylabel([num2str(height),' [mm]'])

figure
plot(xMR_new(14:end-1),yMR_new(14:end-1),'bo')
grid on
hold on
plot(mean(xMR_new(14:end-1)),mean(yMR_new(14:end-1)),'ro')
plot(N/2,M/2,'o','MarkerFaceColor','g')
xlim([0,N])
ylim([0,M])
legend('MR tracked position', 'MR mean position', 'Frame center')
title('MR position inside each RF cineloop')
height = M*pixel_height_RF;
width = N*pixel_width;
% pbaspect([width height 1])
xlabel('[mm]')
ylabel('[mm]')
xticklabels(num2cell((str2double(xticklabels).*pixel_width)))
yticklabels(num2cell((str2double(yticklabels).*pixel_height_RF)))

% %flip 11/05
% R = flip(Reconstruction,2);
% [Mr,Nr] = size(R);
% x_MR = 1+(Nr-xMR_plot(14:end-7));
% y_MR = yMR_plot(14:end-7);
% x_Ve = 1 +(Nr-xMR_plot_ref(14:end-7));
% x_Gr = x_Ve;
% y_Ve = yMR_plot_ref2(14:end-7);
% y_Gr = yMR_plot_ref(14:end-7);
% %convert robot pose
% m = mean(y_img);
% for i=1:length(y_img)
%     y_img(i) = m+sign(m-y_img(i))*abs(m-y_img(i));
% end
% x_pose = x_img./pixel_width;
% x_pose = 1+(Nr-x_pose);
% y_pose = y_img./pixel_height_B;
% y_pose = y_pose+10;
% % x_pose = x_pose-4;
% 
% figure
% imagesc(R);
% colormap('gray');
% hold on
% plot(x_MR,y_MR,'y')
% plot(x_Gr,y_Gr,'g')
% plot(x_Ve,y_Ve,'r')
% plot(x_pose(14:end-6),y_pose(14:end-6),'c')
% height = round(Mr*pixel_height_B);
% width = round(Nr*pixel_width);
% pbaspect([width height 1])
% title('Reconstruction of travelled path. Method 1')
% xlabel([num2str(width),' [mm]'])
% ylabel([num2str(height),' [mm]'])
% legend('Tracked MR','Ground truth','Segmented vessel','Robot pose')
% xticklabels(num2cell(round(str2double(xticklabels).*pixel_width)))
% yticklabels(num2cell(round(str2double(yticklabels).*pixel_height_B)))

%% reconstruction of pose
%to reconstruct pose of joint 5 use direct kinematic functions
%in example the code to reconstruct roll,pitch and yaw of joint 6 is
[M,N] = size(joints);
angles = zeros(M,3);
for i = 1:M
    T = directkinematic_pane(joints(i,:));
    angles(i,:) = tr2rpy_pane(T);
end
angles = angles.*(180/pi);
figure
subplot(3,1,1)
plot(angles(:,1),'b')
grid on
title('Roll angle of joint 6')
ylabel('angle [°]')
ylim([-180,180])
xlabel('sample []')
subplot(3,1,2)
plot(angles(:,2),'b')
grid on
title('Pitch angle of joint 6')
ylabel('angle [°]')
ylim([-90,90])
xlabel('sample []')
subplot(3,1,3)
plot(angles(:,3),'b')
grid on
title('Yaw angle of joint 6')
ylabel('angle [°]')
ylim([-180,180])
xlabel('sample []')
%on joint 6 we can see magnet rotation or vibrations
%to reconstruc US orientation we need to reconstruct joint 5 orientation so
%we need to modify directkinematic_pane and DH_parameters_RV3SB_pane to
%work only on the first 5 joints
%test imposing zero values to joint 6
joints(:,6) = 0; 
angles = zeros(M,3);
for i = 1:M
    T = directkinematic_pane(joints(i,:));
    angles(i,:) = tr2rpy_pane(T);
end
angles = angles.*(180/pi);
figure
subplot(3,1,1)
plot(angles(:,1),'b')
grid on
title('Roll angle of joint 5')
ylabel('angle [°]')
ylim([-180,180])
xlabel('sample []')
subplot(3,1,2)
plot(angles(:,2),'b')
grid on
title('Pitch angle of joint 5')
ylabel('angle [°]')
ylim([-90,90])
xlabel('sample []')
subplot(3,1,3)
plot(angles(:,3),'b')
grid on
title('Yaw angle of joint 5')
ylabel('angle [°]')
ylim([-180,180])
xlabel('sample []')

%we should use the pitch angle of joint 5 to rotate the images in order to
%properly merge them
pitch = angles(indexes,2);
% pitch = -(pitch-pitch(end));
p = find(pitch==max(pitch));
pitch = -(pitch-max(pitch));
pitch(p(end):end) = -pitch(p(end):end);


%Method 2: rf2bmode of merge of cineloop mean frame
[M,N,P] = size(I);
Reconstruction2_rot = zeros(M,N);
start_x = 1;
start_y = 1;
% ref_coord_row = zeros(1,cineloop_nr);
% ref_coord_column = zeros(1,cineloop_nr);
for k = 1:cineloop_nr
    if(k>1)
        start_x = start_x + x_orig(k-1);
        start_y = start_y + y_orig(k-1);
    end
%     ref_coord_column(cineloop_nr-k+1) = start_x;
%     ref_coord_row(cineloop_nr-k+1) = start_y;
    frame = meanI(:,:,cineloop_nr-k+1);
    frame = imrotate(frame,pitch(cineloop_nr-k+1),'crop');
    [M,N] = size(frame);
    if((start_x<=0)||(start_y<=0))
        [Mr,Nr] = size(Reconstruction2_rot);
        Ma = Mr;
        Na = Nr;
        start_x_Re = 1;
        start_y_Re = 1;
        if(start_x<=0)
            Na = Nr-start_x+1;
            start_x_Re = abs(start_x)+1;
            start_x = 1;
        end
        if(start_y<=0)
            Ma = Mr-start_y+1;
            start_y_Re = abs(start_y)+1;
            start_y = 1;
        end
        Aux = Reconstruction2_rot;
        Reconstruction2_rot = zeros(Ma,Na);
        Reconstruction2_rot(start_y:start_y+M-1,start_x:start_x+N-1) = frame;
        Reconstruction2_rot(start_y_Re:start_y_Re+Mr-1,start_x_Re:start_x_Re+Nr-1) = Aux;
%         [Ma,Na] = size(Reconstruction2);
%         Aux = zeros(Ma+abs(start_x)+1,Na+abs(start_y)+1);
%         Aux(abs(start_y)+1:abs(start_y)+Ma,abs(start_x)+1:abs(start_x)+Na) = Reconstruction2;
%         Reconstruction2 = Aux;
%         Reconstruction2(1:M,1:N) = meanI(:,:,cineloop_nr-k+1);
%         if(k<10||k>cineloop_nr-40)
%             Reconstruction2(start_y:start_y+M-1,start_x:start_x+N-1) = zeros(M,N); 
%         end
    else
        Reconstruction2_rot(start_y:start_y+M-1,start_x:start_x+N-1) = frame;
%         if(k<10||k>cineloop_nr-40)
%             Reconstruction2(start_y:start_y+M-1,start_x:start_x+N-1) = zeros(M,N); 
%         end
    end
end
Reconstruction2_rot = rf2bmode(Reconstruction2_rot);

%Method 1: merge rf2bmode images of cineloop mean frame
A = rf2bmode(meanI(:,:,1));
Reconstruction_rot = zeros(size(A));
start_x = 1;
start_y = 1;
start_y2 = 1;
%Exp 1 11_05 comment if(k<10||k>cineloop_nr-40) part
for k = 1:cineloop_nr
    A = rf2bmode(meanI(:,:,cineloop_nr-k+1));
    A = imrotate(A,pitch(cineloop_nr-k+1));
    [M,N] = size(A);
    if(k>1)
        start_x = start_x + x_orig(k-1);
        start_y = start_y + y_orig(k-1);
        start_y2 = start_y2 + floor(y_orig(k-1)/4);
        if(start_y2==0)
            start_y2=1;
        end
    end
%     if(k<10||k>cineloop_nr-40)
%         A = zeros(M,N);
%     end
    if((start_x<=0)||(start_y2<=0))
        [Mr,Nr] = size(Reconstruction_rot);
        Ma = Mr;
        Na = Nr;
        start_x_Re = 1;
        start_y_Re = 1;
        if(start_x<=0)
            Na = Nr-start_x+1;
            start_x_Re = abs(start_x)+1;
            start_x = 1;
        end
        if(start_y2<=0)
            Ma = Mr-start_y2+1;
            start_y_Re = abs(start_y2)+1;
            start_y2 = 1;
        end
        Aux = Reconstruction_rot;
        Reconstruction_rot = zeros(Ma,Na);
        Reconstruction_rot(start_y2:start_y2+M-1,start_x:start_x+N-1) = A;
        Reconstruction_rot(start_y_Re:start_y_Re+Mr-1,start_x_Re:start_x_Re+Nr-1) = Aux;
        
%         [Ma,Na] = size(Reconstruction);
%         Aux = zeros(Ma+abs(start_x)+1,Na+abs(start_y2)+1);
%         Aux(abs(start_y2)+1:abs(start_y2)+Ma,abs(start_x)+1:abs(start_x)+Na) = Reconstruction;
%         Reconstruction = Aux;
%         Reconstruction(1:M,1:N) = A;
    else
        Reconstruction_rot(start_y2:start_y2+M-1,start_x:start_x+N-1) = A;
    end
end

figure
imagesc(Reconstruction_rot);
colormap('gray');
[Mr,Nr] = size(Reconstruction_rot);
height = round(Mr*pixel_height_B);
width = round(Nr*pixel_width);
pbaspect([width height 1])
title('Reconstruction of travelled path. Method 1')
xlabel([num2str(width),' [mm]'])
ylabel([num2str(height),' [mm]'])

figure
imagesc(Reconstruction2_rot);
colormap('gray');
[Mr,Nr] = size(Reconstruction2_rot);
height = round(Mr*pixel_height_B);
width = round(Nr*pixel_width);
pbaspect([width height 1])
title('Reconstruction of travelled path. Method 2')
xlabel([num2str(width),' [mm]'])
ylabel([num2str(height),' [mm]'])

%% vessel segmentation
% [Mr,Nr] = size(Reconstruction);
% filter = [-1;1];
% Rec_filtered = abs(imfilter(Reconstruction,filter));
% filter2 = ones(5,5)./25;
% % Rec_filtered = abs(imfilter(Rec_filtered,filter2));
% %we should now detect the first peak at th value, without considering the
% %possible subsequent peaks
% th = 15;
% % Rec_filtered = abs(imfilter(Rec_filtered>th,filter));
% Rec_filtered = double(Rec_filtered);%>th);
A = rf2bmode(meanI(:,:,1));
[M,N] = size(A);
A(1:3,:) = 0;
A(M-2:M,:) = 0;
filter = [zeros(1,60); ones(1,60); zeros(1,60)];
A_fil = imfilter(A,filter);
[pks,locs] = findpeaks(A_fil(:,10));
x = 10;
y = locs(pks==max(pks));
pks(pks==max(pks)) = 0;
y = (y+locs(pks==max(pks)))/2;
figure
plot(pks)
figure
plot(locs)
figure
imagesc(A_fil)
colormap('gray')
% Rec_filtered = sobel(A,th);
% Rec_filtered(1:10,:) = 0;
% Rec_filtered(M-5:M,:) = 0;
% Rec_filtered(:,1:3) = 0;
% Rec_filtered(:,N-3:N) = 0;
% figure
% % plot(Rec_filtered(:,5));
% imshow(Rec_filtered);
% Rec_filtered = double(Rec_filtered>0);
% for i=1:N
%     [pks,locs] = findpeaks(Rec_filtered(:,i));
%     if(length(locs)>1)
%         Rec_filtered(locs(2):M,i) = 0;
%     end
% %     start = 1;
% %     for j=2:Mr
% %         if(Rec_filtered(j-1,i)~=Rec_filtered(j,i))
% %             start = j;
% %             break
% %         end
% %     end
% %     Rec_filtered(start:Mr,i) = 0;
% end
% % Rec_filtered(round(Mr/2):Mr,:) = 0;
% [vessel_y,vessel_x] = find(Rec_filtered==1);
start_x = 1;
start_y = 1;
for i = 1:length(ref_x_vessel)
    m = ref_x_vessel(i);
    n = ref_y_vessel(i);
    vessel_x(start_x:start_x+m-1) = vessel_x(start_x:start_x+m-1)+ref_coord_column(cineloop_nr-i+1)*ones(m,1);
    vessel_y(start_y:start_y+n-1) = vessel_y(start_y:start_y+n-1)+ref_coord_row(cineloop_nr-i+1)*ones(n,1);
    start_x = start_x+m;
    start_y = start_y+n;
end
figure
plot(vessel_x,vessel_y)
[Mr,Nr] = size(Reconstruction2);
xlim([0 Nr])
ylim([0 Mr])
ax = gca;
ax.YDir = 'reverse';

% %% error
% %exp1 11/05
% x_pose = 1+(Nr-x_pose);
% error_x = abs(xMR_plot(14:end)-xMR_plot_ref(14:end)).*pixel_width;
% error_y = abs(yMR_plot(14:end)-yMR_plot_ref(14:end)).*pixel_height_B;
% error_x_off = abs(offline_xMR'-xMR_plot_ref(14:end)).*pixel_width;
% error_y_off = abs(offline_yMR-yMR_plot_ref(14:end)).*pixel_height_B;
% MR_diam = 0.55;
% th_x = (MR_diam/2).*ones(1,length(error_x));
% th_y = (MR_diam/2).*ones(1,length(error_y));
% errorx = error_x(end-24:end);
% errory = error_y(end-24:end);
% error_x = error_x(1:end-24);
% error_y = error_y(1:end-24);
% error_x2 = abs(x_pose(14:end-7)-xMR_plot_ref(14:end-7)).*pixel_width;
% error_y2 = abs(y_pose(14:end-7)'-yMR_plot_ref(14:end-7)).*pixel_height_B;
% n = length(error_x);
% sample = n + (1:length(errorx))-1;
% sample1 = 1:length(error_x);
% figure
% subplot(2,1,1)
% plot(sample1,error_x,'b')
% hold on
% grid on
% plot(sample,errorx,'r')
% plot(th_x,'k--')
% xlim([0 70])
% ylim([0 3])
% ylabel('[mm]')
% title('MR tracking error on x')
% subplot(2,1,2)
% plot(sample1,error_y,'b')
% hold on
% grid on
% plot(sample,errory,'r')
% plot(th_y,'k--')
% xlim([0 70])
% ylim([0 3])
% ylabel('[mm]')
% title('MR tracking error on y')
% legend('inside image plane','outside image plane','one radius')
% 
% figure
% subplot(1,2,1)
% boxplot(error_x','Labels',{['Online: µ = ',num2str(mean(error_x)), ' mm']})
% grid on
% ylabel('[mm]')
% ylim([0 3])
% title('Error on x')
% subplot(1,2,2)
% boxplot(error_y','Labels',{['Online: µ = ',num2str(mean(error_y)), ' mm']})
% grid on
% ylabel('[mm]')
% ylim([0 3])
% title('Error on y')
% 
% E = [error_x,error_y'];
% figure
% boxplot(E,'Labels',{['Error on x: µ = ',num2str(mean(error_x)), ' mm'],['Error on y: µ = ',num2str(mean(error_y)), ' mm']})
% grid on
% ylabel('[mm]')
% ylim([0 2])
% title('Tracking error')
% 
% error_x_off = error_x_off(1:end-24);
% error_y_off = error_y_off(1:end-24);
% E = [error_x_off,error_y_off'];
% figure
% boxplot(E,'Labels',{['Error on x: µ = ',num2str(mean(error_x_off)), ' mm'],['Error on y: µ = ',num2str(mean(error_y_off)), ' mm']})
% grid on
% ylabel('[mm]')
% ylim([0 2])
% title('Tracking error offline')
% 
% 
% figure
% subplot(2,1,1)
% plot(error_x2,'b')
% hold on
% grid on
% plot(th_x,'k--')
% xlim([0 70])
% ylim([0 3])
% ylabel('[mm]')
% title('MR tracking error on x')
% subplot(2,1,2)
% plot(error_y2,'b')
% hold on
% grid on
% plot(th_y,'k--')
% xlim([0 70])
% ylim([0 3])
% ylabel('[mm]')
% title('MR tracking error on y')
% legend('online','one radius')
% figure
% subplot(1,2,1)
% boxplot(error_x2','Labels',{['Online: µ = ',num2str(mean(error_x2)), ' mm']})
% grid on
% ylabel('[mm]')
% ylim([0 3])
% title('Error on x')
% subplot(1,2,2)
% boxplot(error_y2','Labels',{['Online: µ = ',num2str(mean(error_y2)), ' mm']})
% grid on
% ylabel('[mm]')
% ylim([0 3])
% title('Error on y')


% %exp1 30/06
% error_x = abs(xMR_plot(25:end)-xMR_plot_ref(25:end)).*pixel_width;
% error_y = abs(yMR_plot(25:end)-yMR_plot_ref(25:end)).*pixel_height_B;
% MR_diam = 1.25;
% th_x = (MR_diam/2).*ones(1,length(error_x));
% th_y = (MR_diam/2).*ones(1,length(error_y));
% figure
% subplot(2,1,1)
% plot(error_x,'b')
% hold on
% grid on
% plot(th_x,'k--')
% xlim([0 50])
% ylim([0 2])
% ylabel('[mm]')
% title('MR tracking error on x')
% subplot(2,1,2)
% plot(error_y,'b')
% hold on
% grid on
% plot(th_y,'k--')
% xlim([0 50])
% ylim([0 2])
% ylabel('[mm]')
% title('MR tracking error on y')
% legend('online','one radius')
% 
% figure
% subplot(1,2,1)
% boxplot(error_x','Labels',{['Online: µ = ',num2str(mean(error_x)), ' mm']})
% grid on
% ylabel('[mm]')
% ylim([0 2])
% title('Error on x')
% subplot(1,2,2)
% boxplot(error_y','Labels',{['Online: µ = ',num2str(mean(error_y)), ' mm']})
% grid on
% ylabel('[mm]')
% ylim([0 2])
% title('Error on y')
% 
% save('Error_exp1_bignoflux_30_06.mat','error_x','error_y')

% %exp2 30/06
% error_x = abs(xMR_plot(25:end-10)-xMR_plot_ref(25:end-10)).*pixel_width;
% error_y = abs(yMR_plot(25:end-10)-yMR_plot_ref(25:end-10)).*pixel_height_B;
% MR_diam = 1.25;
% th_x = (MR_diam/2).*ones(1,length(error_x));
% th_y = (MR_diam/2).*ones(1,length(error_y));
% figure
% subplot(2,1,1)
% plot(error_x,'b')
% hold on
% grid on
% plot(th_x,'k--')
% xlim([0 65])
% ylim([0 2])
% ylabel('[mm]')
% title('MR tracking error on x')
% subplot(2,1,2)
% plot(error_y,'b')
% hold on
% grid on
% plot(th_y,'k--')
% xlim([0 65])
% ylim([0 2.5])
% ylabel('[mm]')
% title('MR tracking error on y')
% legend('online','one radius')
% 
% figure
% subplot(1,2,1)
% boxplot(error_x','Labels',{['Online: µ = ',num2str(mean(error_x)), ' mm']})
% grid on
% ylabel('[mm]')
% ylim([0 2])
% title('Error on x')
% subplot(1,2,2)
% boxplot(error_y','Labels',{['Online: µ = ',num2str(mean(error_y)), ' mm']})
% grid on
% ylabel('[mm]')
% ylim([0 2])
% title('Error on y')

% save('Error_exp2_bignoflux_30_06.mat','error_x','error_y')

% %exp3 30/06
% error_x = abs(xMR_plot(106:end)-xMR_plot_ref(106:end)).*pixel_width;
% error_y = abs(yMR_plot(106:end)-yMR_plot_ref(106:end)).*pixel_height_B;
% MR_diam = 1.25;
% th_x = (MR_diam/2).*ones(1,length(error_x));
% th_y = (MR_diam/2).*ones(1,length(error_y));
% figure
% subplot(2,1,1)
% plot(error_x,'b')
% hold on
% grid on
% plot(th_x,'k--')
% xlim([0 45])
% ylim([0 2])
% ylabel('[mm]')
% title('MR tracking error on x')
% subplot(2,1,2)
% plot(error_y,'b')
% hold on
% grid on
% plot(th_y,'k--')
% xlim([0 45])
% ylim([0 2.5])
% ylabel('[mm]')
% title('MR tracking error on y')
% legend('online','one radius')
% 
% figure
% subplot(1,2,1)
% boxplot(error_x','Labels',{['Online: µ = ',num2str(mean(error_x)), ' mm']})
% grid on
% ylabel('[mm]')
% ylim([0 2])
% title('Error on x')
% subplot(1,2,2)
% boxplot(error_y','Labels',{['Online: µ = ',num2str(mean(error_y)), ' mm']})
% grid on
% ylabel('[mm]')
% ylim([0 2])
% title('Error on y')
% 
% save('Error_exp3_bigflux_30_06.mat','error_x','error_y')

%exp4 30/06
error_x = abs(xMR_plot(10:end)-xMR_plot_ref(10:end)).*pixel_width;
error_y = abs(yMR_plot(10:end)-yMR_plot_ref(10:end)).*pixel_height_B;
MR_diam = 0.55;
th_x = (MR_diam/2).*ones(1,length(error_x));
th_y = (MR_diam/2).*ones(1,length(error_y));
figure
subplot(2,1,1)
plot(error_x,'b')
hold on
grid on
plot(th_x,'k--')
xlim([0 100])
ylim([0 2])
ylabel('[mm]')
title('MR tracking error on x')
subplot(2,1,2)
plot(error_y,'b')
hold on
grid on
plot(th_y,'k--')
xlim([0 100])
ylim([0 2])
ylabel('[mm]')
title('MR tracking error on y')
legend('online','one radius')

figure
subplot(1,2,1)
boxplot(error_x','Labels',{['Online: µ = ',num2str(mean(error_x)), ' mm']})
grid on
ylabel('[mm]')
ylim([0 2])
title('Error on x')
subplot(1,2,2)
boxplot(error_y','Labels',{['Online: µ = ',num2str(mean(error_y)), ' mm']})
grid on
ylabel('[mm]')
ylim([0 2])
title('Error on y')
% 
% save('Error_exp4_smallflux_30_06.mat','error_x','error_y')

% %exp5 30/06
% error_x = abs(xMR_plot(28:end)-xMR_plot_ref(28:end)).*pixel_width;
% error_y = abs(yMR_plot(28:end)-yMR_plot_ref(28:end)).*pixel_height_B;
% MR_diam = 0.55;
% th_x = (MR_diam/2).*ones(1,length(error_x));
% th_y = (MR_diam/2).*ones(1,length(error_y));
% figure
% subplot(2,1,1)
% plot(error_x,'b')
% hold on
% grid on
% plot(th_x,'k--')
% xlim([0 30])
% ylim([0 2])
% ylabel('[mm]')
% title('MR tracking error on x')
% subplot(2,1,2)
% plot(error_y,'b')
% hold on
% grid on
% plot(th_y,'k--')
% xlim([0 30])
% ylim([0 2])
% ylabel('[mm]')
% title('MR tracking error on y')
% legend('online','one radius')
% 
% figure
% subplot(1,2,1)
% boxplot(error_x','Labels',{['Online: µ = ',num2str(mean(error_x)), ' mm']})
% grid on
% ylabel('[mm]')
% ylim([0 2])
% title('Error on x')
% subplot(1,2,2)
% boxplot(error_y','Labels',{['Online: µ = ',num2str(mean(error_y)), ' mm']})
% grid on
% ylabel('[mm]')
% ylim([0 2])
% title('Error on y')
% 
% save('Error_exp5_smallnoflux_30_06.mat','error_x','error_y')

% %exp6 30/06
% error_x = abs(xMR_plot(101:end-5)-xMR_plot_ref(101:end-5)).*pixel_width;
% error_y = abs(yMR_plot(101:end-5)-yMR_plot_ref(101:end-5)).*pixel_height_B;
% MR_diam = 0.55;
% th_x = (MR_diam/2).*ones(1,length(error_x));
% th_y = (MR_diam/2).*ones(1,length(error_y));
% figure
% subplot(2,1,1)
% plot(error_x,'b')
% hold on
% grid on
% plot(th_x,'k--')
% xlim([0 85])
% ylim([0 2.5])
% ylabel('[mm]')
% title('MR tracking error on x')
% subplot(2,1,2)
% plot(error_y,'b')
% hold on
% grid on
% plot(th_y,'k--')
% xlim([0 85])
% ylim([0 2.5])
% ylabel('[mm]')
% title('MR tracking error on y')
% legend('online','one radius')
% 
% figure
% subplot(1,2,1)
% boxplot(error_x','Labels',{['Online: µ = ',num2str(mean(error_x)), ' mm']})
% grid on
% ylabel('[mm]')
% ylim([0 2.5])
% title('Error on x')
% subplot(1,2,2)
% boxplot(error_y','Labels',{['Online: µ = ',num2str(mean(error_y)), ' mm']})
% grid on
% ylabel('[mm]')
% ylim([0 2.5])
% title('Error on y')
% 
% save('Error_exp6_smallnoflux_30_06.mat','error_x','error_y')