qu  clear all
close all
clc

% Get info (including filename) for each file in current directory
list = dir('*.JPG');
number_of_files = size(list);
 
% Define rectangle near centre of cap. Determined EMPIRICALLY.
cap_x = [165, 185, 185, 165];
cap_y = [20, 20, 35, 35];
cap_thres = 0.04; % Any hue value less than a hue of 0.04 as red in these images.

% Define rectangle for overfill case. Determined EMPIRICALLY.
over_x = [165, 185, 185, 165];
over_y = [95, 95, 110, 110];
over_thres = 0.5; % Any V value in HSV colour model less than 0.5 is black in our images.

% Define rectangle for underfill case. Determined EMPIRICALLY.
under_x = [155, 175, 175, 155];
under_y = [160, 160, 175, 175];
under_thres = 0.5; % Any V value in HSV colour model less than 0.5 is black in our images.

% Go through each image file one at a time and test for faults.
for i= 1: number_of_files(1,1)
    
    % Read file and display image
    filename = list(i).name;
    im = imread(filename);
    imshow(im);
    
    % Plot rectangles corresponding to cap, overfill region and underfill
    % region.
    imrect(gca, [cap_x(1), cap_y(1), cap_x(2)-cap_x(1)+1, cap_y(3)-cap_y(1)+1]);
    imrect(gca, [over_x(1), over_y(1), over_x(2)-over_x(1)+1, over_y(3)-over_y(1)+1]);
    imrect(gca, [under_x(1), under_y(1), under_x(2)-under_x(1)+1, under_y(3)-under_y(1)+1]);
    
    % Convert image from RGB to HSV for better definition of hues (i.e.,
    % colours such as red.
    % Also allows separation of value component (brightness) which can be
    % tested for "blackness".
    hsv = rgb2hsv(im);
    hue = hsv(:,:,1);
    value = hsv(:,:,3);
    
    % Check if cap present or not. If cap present, then cap region should
    % be red, which is a hue value close to zero (<cap_thres for our images).
    capped = 'Capped. '; % Assume capped by default
    cap_hue = mean2( hue( cap_y(1):cap_y(3), cap_x(1):cap_x(2)) ); % Compute mean hue in cap rectangle
    if (cap_hue > cap_thres)
        capped = 'Not capped.';
    end

    filling = 'Normally filled.'; % Assume normal filling by default
    % For overfilling, it is best to process the value component.
    % Black is close to zero or <over_thres for our images.
    %
    %**********************************************************************
    % WRITE CODE HERE THAT (1) computes the mean V value in rectangle
    % defined by over_x and over_y, (2) checks if mean value is less than
    % over_thres, and (3) sets message accordingly if mean is less than
    % over_thres (i.e., filling = 'Overfilled.')
    %**********************************************************************
    
    V = mean2( value( over_y(1):over_y(3), over_x(1):over_x(2)) );
     if (V < over_thres)
        filling = 'Overfilled.';
    end
    
    % For underfilling, it is best to process the value component.
    % Black is close to zero or <under_thres for our images.
    %
    %**********************************************************************
    % WRITE CODE HERE THAT (1) only executes steps (2) - (4) if bottle is
    % not overfilled, (2) computes the mean V value in rectangle defined
    % by under_x and under_y, (3) checks if mean value is greater than
    % under_thres, and (3) sets message accordingly if mean is greater than
    % under_thres (i.e., filling = 'Underfilled.')
    % For step 1, you may wish to use the MATLAB function strcmp to check
    % if filling is equal to 'Overfilled.'
    %**********************************************************************
    V1 = mean2( value( under_y(1):under_y(3), under_x(1):under_x(2)) );
     if (V1 > under_thres)
        filling = 'Underfilled.';
    end
    
    message = strcat(capped, filling);
    title(message)
    xlabel('Press mouse button to continue')
    % Compute mean value in ROI defined by rectangle
    waitforbuttonpress
end


