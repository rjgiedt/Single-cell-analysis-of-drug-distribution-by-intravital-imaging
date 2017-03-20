%% Initial setup

clear all; clc; close all;

%% Input Video (Cell location)
input_video_fname = 'Input.avi';
              
%Prompt for input video file name (use above name as default)
[input_video_fname,input_video_path] = uigetfile('*.avi',...
    'Select the cell/nuclear boundary input video',input_video_fname);
if isequal(input_video_fname,0) || isequal(input_video_path,0)
   disp('User selected Cancel')
   return;
else
   fprintf('User selected "%s" as input video\n',input_video_fname);
end

%% Input Video (Cellular drug Concentrations)
input_video_fname2 = 'Input2.avi';
              
%Prompt for input video file name (use above name as default)
[input_video_fname2,input_video_path2] = uigetfile('*.avi',...
    'Select the drug input video',input_video_fname2);
if isequal(input_video_fname2,0) || isequal(input_video_path2,0)
   disp('User selected Cancel')
   return;
else
   fprintf('User selected "%s" as input video\n',input_video_fname2);
end

%% Prompt for filtered video output name
output_video_fname = 'Output.avi';
[output_video_fname,output_video_path] = uiputfile('*.avi','Input output video file name');

if isequal(output_video_fname,0) || isequal(output_video_path,0)
   disp('User selected Cancel')
   return;
else
   fprintf('User selected "%s" as output video file\n',output_video_fname);
end

%% Read in selected videos 
movie=VideoReader(input_video_fname);  %Read in cell boundary video
nframes = movie.NumberOfFrames;  %Determine the total number of frames

% Read in video of drug diffusion 
movie2=VideoReader(input_video_fname2); %Read in Drug diffusion video
nframes2 = movie2.NumberOfFrames; %Determine total number of frames for drug diffusion video

%% Display Error for Inconsistent frame numbers
if nframes~=nframes2
    warning('Videos must be equal Length... video length mismatch')
end
    
%% Define Image and allow user to define correct gamma value based on image output/ filtering

goAgain = true;
while goAgain

%close previous windows
close all;

%Query user for thresholding values
prompt={'Enter gamma value to adjust image', 'Enter video frame to check',...
    'Enter Iteration Number', 'Enter Power Value',...
    'Enter Filtering Disk Size (px)', 'Enter Object Minimum (px)'};

title='Image thresholding calculator'; 

%Define variables
answer=inputdlg(prompt,title);
gamma = str2num(answer{1});
videoframenumber = str2num(answer{2});
iteration = str2num(answer{3});
power = str2num(answer{4}); 
disk_size = str2num(answer{5});
cutoff = str2num(answer{6}); 

%Import image and filter using values described 
I = read(movie,videoframenumber);
Threshold_Image = I; %rgb2gray(
Region_Image_Adj = imadjust(Threshold_Image);
negative = imadjust(Threshold_Image, [0 1], [1 0]);
gamma_adjust = imadjust(negative, [], [], gamma);

%Run and display thresholding sequence
[B,T] = Minimax_output(gamma_adjust, iteration, power);

% Filter and display final output of original image

se = strel('disk',disk_size); % Set disk size

Io = imopen(B,se);
Ie = imerode(B, se);
Iobr = imreconstruct(Ie, B);
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
Iobrcbr = im2double(Iobrcbr);
Iobrcbr = imadjust(Iobrcbr, [0,1], [1,0]);
labeled_filled = imfill(Iobrcbr, 'holes');

%Conduct Size Filtering dependent on user input
size_filtering = bwlabel(labeled_filled);
s = regionprops(size_filtering, 'Area');
area_values = [s.Area];
idx = find((cutoff <= area_values) & (200000 >= area_values));

final_image = ismember(size_filtering, idx);

labeled = bwlabel(final_image);

% Use bwlabel to color code perimeters of analyzed regions 
doubled = im2double(labeled);
negative2 = imadjust(doubled, [0,1], [1,0]);
labeled = bwlabel(negative2);
color_labeled = label2rgb(labeled, 'hsv', 'w', 'shuffle');

% Label Objects and display overlay to see relative region selections
bw4_perim = bwperim(negative2);
overlay1 = imoverlay(gamma_adjust, bw4_perim, [.3 1 .3]);
overlay2 = imoverlay(I, bw4_perim, [.3 1 .3]);

%% Create Subplot of output for user to evaluate
figure; 
imshow(I)
%title('Original Image');

figure;
imshow(gamma_adjust)
%title('Gamma Adjusted Image');

figure;
imshow(color_labeled)
%title('Color Labeled Image');

figure;
imshow(overlay1)
%title('Negative Image Overlay');

figure;
imshow(overlay2)
%title('Original Image Overlay');


%% Prompt user if they would like to repeat process with a new gamma?  
prompt={'Would you like to retry a new value for gamma (Y/N)?'};
 title = 'Repeat Analysis?';   
 answer = inputdlg(prompt,title);
 res = answer{1};
 goAgain = isequal(upper(res),'Y');
end

%Query user for output video framerate
prompt2={'Enter Frame Rate of output video'};
title='Analyzed Video Output'; 
answer = inputdlg(prompt2,title);
res = answer{1};
frame_rate = str2num(answer{1}); 

%% Create a filtered video from the above defined conditions
close all;

% Create wait bar for thresholding
h = waitbar(0,'Frames Completed','Name','Thresholding Images',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)

%Construct a VideoWriter object and view its properties. 
writerObj = VideoWriter(output_video_fname);

%Set the frame rate to 20 frames per second.
writerObj.FrameRate = frame_rate;  
open(writerObj);

%Run loop to filter all video frames 
for k = 1:nframes
   P = read(movie,k);
   Q = rgb2gray(read(movie2,k));
 
close all
 
Threshold_Image = P;%rgb2gray(
Region_Image_Adj = imadjust(Threshold_Image);
negative = imadjust(Threshold_Image, [0 1], [1 0]);
gamma_adjust = imadjust(negative, [], [], gamma);

%Run and display thresholding sequence
[B,T] = Minimax_output2(gamma_adjust, iteration, power);

% Filter and display final output of original image

se = strel('disk',disk_size); % Set disk size

Io = imopen(B,se);
Ie = imerode(B, se);
Iobr = imreconstruct(Ie, B);
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
Iobrcbr = im2double(Iobrcbr);
Iobrcbr = imadjust(Iobrcbr, [0,1], [1,0]);
labeled_filled = imfill(Iobrcbr, 'holes');

%Conduct Size Filtering dependent on user input
size_filtering = bwlabel(labeled_filled);
s = regionprops(size_filtering, 'Area');
area_values = [s.Area];
idx = find((cutoff <= area_values) & (200000 >= area_values));

final_image = ismember(size_filtering, idx);

[labeled num] = bwlabel(final_image);
number(k,1) = num;

% Generate statistics for pk analysis/ tracking uploading
stats = regionprops(labeled, 'centroid'); % Centroid generation
stats2 = regionprops(labeled, Threshold_Image, 'MeanIntensity'); %Intensity of original images 
stats3 = regionprops(labeled, Q, 'MeanIntensity'); %Obtain Drug Intensity information here
stats4 = regionprops(labeled, 'Area'); % Track Area

% Calculate X and Y Centroid locations for each particle and store in their
% own cells for graphing

    for i=1:num
    locstatsX{i,k} = (stats(i).Centroid(1));
    locstatsY{i,k} = (stats(i).Centroid(2));
    locstatsAmp{i,k} = (stats2(i).MeanIntensity);
    drugstats{i,k} = stats3(i).MeanIntensity;
    area_stats{i,k} = stats4(i).Area;
    end
    
    
% Assemble Information for Danuser program
    for i=1:num
    locstatsX2{i,k} = (stats(i).Centroid(1));
    locstatsY2{i,k} = (stats(i).Centroid(2));
    locstatsAmp2{i,k} = (stats3(i).MeanIntensity);
    end
    
% Use bwlabel to color code perimeters of analyzed regions
doubled = im2double(labeled);
negative2 = imadjust(doubled, [0,1], [1,0]);
labeled = bwlabel(negative2);

% Generate Overlay Videos for analysis
bw4_perim = bwperim(negative2);
overlay1 = imoverlay(gamma_adjust, bw4_perim, [.3 1 .3]);
overlay2 = imoverlay(I, bw4_perim, [.3 1 .3]);

% Write Cell Borders Output Video    
 inputframe = im2double(overlay1);
 frame = inputframe;
 writeVideo(writerObj,frame);  
    
% Alert user to status of loop
 % Check for Cancel button press
    if getappdata(h,'canceling')
        break
    end
    % Report current estimate in the waitbar's message field
    waitbar(k/nframes)    
end
close(writerObj);

%% Close all windows including status bar 
set(0,'ShowHiddenHandles','on');
delete(get(0,'Children'));
set(0,'ShowHiddenHandles', 'off');

%% Assemble Information for External Linking Program

for i = 1:nframes
    %Structures logging centroid data
    movieInfo(i,1).xCoord(:,1) = cell2mat(locstatsX2(:,i));
    movieInfo(i,1).yCoord = cell2mat(locstatsY2(:,i));
    movieInfo(i,1).amp = cell2mat(locstatsAmp2(:,i));
   
    %Structures logging Stdevs (set to zero)
    movieInfo(i,1).xCoord(:,2) = zeros(size(cell2mat(locstatsX2(:,i))));
    movieInfo(i,1).yCoord(:,2) = zeros(size(cell2mat(locstatsY2(:,i))));
    movieInfo(i,1).amp(:,2) = zeros(size(cell2mat(locstatsAmp2(:,i))));
    %movieInfo(i,1).drug(:,2) = zeros(size(cell2mat(drugstats(:,i))));
end


%% Calculate and graph intracellular drug concentrations

%Initial Graph (no filtering of regions)
%numbers = cell2mat(drugstats)
for i=1:nframes
Frame_Averages(i)= mean(cell2mat(drugstats(:,i)));   
Total_Area(i) = sum(cell2mat(area_stats(:,i)));
standard_deviation(i) = mean(std2(cell2mat(drugstats(:,i))));
end
figure; 
plot(1:nframes, Frame_Averages, 1:nframes, standard_deviation);
%title('Average Curves');

%% Plot number of cells making some minimum drug concentration vs. time

prompt3={'Enter Value for minimum therapuetic drug concentration'};
title='Minimum Therapuetic Dosing'; 
answer = inputdlg(prompt3,title);
res = answer{1};
thera_minimum = str2num(answer{1}); 

therapuetic_value = thera_minimum;
 
for i=1:nframes
     for j = 1:number(i)
         average(i,j) = mean(movieInfo(i,1).amp(j,1));
     end
 average2(i) = mean(average(i,:));      
 end

for i = 1:nframes
        dose(i) = numel(average(i,:), average(i,:) > therapuetic_value);
        total_count(i) = numel(average(i,:));
        percentage(i) = dose(i)/total_count(i);
end

figure; plot(percentage)

%% Plot min-max curves

for i = 1:nframes
    min_value(i) = min(cell2mat(drugstats(:,i)));
    max_value(i) = max(cell2mat(drugstats(:,i)));
end

figure;
plot(1:nframes, min_value, 1:nframes, max_value);
%title('Min-Max Curves');


%% Save Output from object recognition phase in mfile
save('drug_movie_data.mat');
