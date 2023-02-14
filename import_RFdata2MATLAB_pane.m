function [RF_DATA, HEADER, transducer_code,flag]=import_RFdata2MATLAB_pane(filename)
%% Function which allows to import RF data and other data acquisition parameters needed for imaging 

%% IN:
% filename - *.bin binary file recorded using ArtUS device and artus_rf_test sample
% DIR - directory of the file

%% OUT:
% RF_DATA - cell format (1 x K - number of RF frames and inside N (number of RF samples) X M (number of RF rows)) 
% HEADER - RF data acquisition parameters: 
% HEADER{1,1}.tx_frequency - ultrasonic wave (transmission) frequency 
% HEADER{1,1}.Length_of_RF_row - number of samples in each RF row
% HEADER{1,1}.Number_of_RF_rows - number of RF scanning lines in a window
% HEADER{1,1}.Sampling_period_ns - sampling period in ns
% HEADER{1,1}.BIT_ADC - number of ADC bits
% HEADER{1,1}.beam_x - start point x coordinates for each ultrasound beam in cm 
% HEADER{1,1}.beam_y - start point y coordinates for each ultrasound beam in cm
% HEADER{1,1}.angle - angle in radians (the angle is given relative to  the perpendicular to the center of the probe’s surface)

% transducer_code - type, bandwidth limits (HF...LF), width of array, code of manufacturer
% flag - indicates sucessfull file reading

fid=fopen([filename]);                                                       

if (fid~=-1)                                                               %% Error check if file is empty 
 File_size=length(fread(fid));                                             %% Estimates the size of file in bytes 
 fclose(fid);
 transducer_code=filename(21:end-4);                                       %% Char, transducer code: type, bandwidth limits (HF...LF), width of array, code of manufacturer
 flag=1;
else  
 errordlg('Selected empty file')
 flag=0;
 RF_DATA=[]; 
 HEADER=[]; 
 transducer_code=[];
 return;
end

fid=fopen([filename]);
number_of_header_els=7;                                                    %% Number of header parameters for each frame 

%% loop control variables 
size_in_bytes_global=0;                                                      
stop=true;
i=1;
%% --------------------------------------------------------------------------------------------------------------------------------------------

%% --------------------------------------------------------------------------------------------------------------------------------------------
%% Loop for RF frames reading and reshaping and construction of header for each frame (in case of scanning parameter changes, i.e. RF window size) 
while (stop)
   
A = fread(fid,number_of_header_els,'int32');                               %% Reads *.bin file recorded using ArtUS device and artus_rf_test sample 

header.tx_frequency=A(1);                                                  %% Transmission frequency, in Hz
header.frame_rate=A(2)/100;                                                %% The actual number of frames per second, since streaming started 
header.Length_of_RF_row=A(3);                                              %% Number of RF signal samples (1 symbol in data file)
header.Number_of_RF_rows=A(4);                                             %% Number of RF signal rows in a predefined RF Window
header.Sampling_period_ns=A(5);                                            %% Sampling_period, nanoseconds
header.BIT_ADC=A(6);                                                       %% No of BIT's per sample (ADC digitization)
start_depth=A(7);                                                          %% Scanning start depth in mm (offset for each beam)

%% ErrorCheck for binary file reading mistakes
if ((header.BIT_ADC~=16)&&(header.BIT_ADC~=32))
     errordlg('Wrong structure of the file')
     flag=0;
     RF_DATA=[]; 
     HEADER=[]; 
     transducer_code=[];
     return; 
else
     flag=1;
end    

%% Start position and orientation of each beam
Start_point_position_and_orient=fread(fid,3*header.Number_of_RF_rows,'int32'); 

%% Relative angle for each beam (0 for linear array)
header.angle=Start_point_position_and_orient(3:3:end)./1000000;             %% Specifies the angle of the ultrasonic beam’s direction. The angle is given relative to  the perpendicular to the center of the probe’s surface.
%% Start position (x,y) of each beam
header.beam_x=Start_point_position_and_orient(1:3:end)/10000+(start_depth/10*sin(header.angle));  %% Start x coordinates of each beam, cm
header.beam_y=Start_point_position_and_orient(2:3:end)/10000+(start_depth/10*cos(header.angle));  %% Start y coordinates of each beam, cm

header.Speed_of_sound=1540;                                                %% Speed of Sound, m/s

HEADER{i}=header;

size_of_frame=header.Length_of_RF_row*header.Number_of_RF_rows;            %% Number of elements of single RF frame  
frame_to_reshape=((fread(fid,size_of_frame,'int16')));                     %% Read single RF frame vector (int16 data type)

%% RF stream reshaped to matrix (all the collected full frames)
RF_DATA{i}=reshape(frame_to_reshape,header.Length_of_RF_row,header.Number_of_RF_rows);  

size_in_bytes_run=(number_of_header_els*4)+(3*header.Number_of_RF_rows*4)+(size_of_frame*2);
size_in_bytes_global=size_in_bytes_global+size_in_bytes_run;
i=i+1;
if (File_size-size_in_bytes_global)>=size_in_bytes_run
  stop=true;
else
  stop=false;
end
clear frame_to_reshape 
end





