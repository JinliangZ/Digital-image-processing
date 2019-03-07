function [im2, a] = autolevel_16( fname )
% autolevel_16
% automatically levels background of an image.
%
%
% USAGE: [im2, a] = autolevel_16( fname ) 
% where 'fname' is a string representing
% the filename of the file containing the image
% 'im2' is the levelled image
% 'a' is a vector of the 6 unknown coefficients.

% Read in image and display it
im = imread( fname );

% block size: bm*bn
bm = 16;
bn = 16;
[rows, cols] = size(im); 
n_points = (rows/bm) * (cols/bn);

x = zeros(1, n_points) ;
y = zeros(1, n_points) ;
I = zeros(1, n_points) ;
idx = 1;

% Automatically select points 
for di = 1 : rows/bm
    for dj = 1 : cols/bn 
        imblock = im( (di-1)*bm+1 : di*bm , (dj-1)*bn+1 : dj*bn ) ;
 
        % find the minimum point's coordination
        [~, tI] = min(imblock(:)) ;
        [ty, tx] = ind2sub(size(imblock), tI) ;
        tv = imblock(ty, tx) ;

        x(idx) = tx + (dj-1)*bn;
        y(idx) = ty + (di-1)*bm;
        I(idx) = tv;
        
        idx = idx+1;
    end
end
 

% Plot points  
figure;
subplot(1,2,1);
imshow(im);
hold on ;
plot(x, y, 'y+') ;
hold off;

% Fit data at selected points to background function
%  Solve least-squares problem: [C]{a} = {k} using the
%    \ operator, i.e., {a} = [C]\{k}
%  First, compute elements of the matrix [C]
N = length(x);
Sx = sum(x);
Sy = sum(y);
Sxx = sum(x.*x);
Syy = sum(y.*y);
Sxy = sum(x.*y);
Sxxx = sum(x.^3);
Sxxy = sum(x.*x.*y);
Sxyy = sum(x.*y.*y);
Syyy = sum(y.^3);
Sxxxx = sum(x.^4);
Sxxxy = sum(y.*x.^3);
Sxxyy = sum(x.*x.*y.*y);
Sxyyy = sum(x.*y.^3);
Syyyy = sum(y.^4);

C = [N    Sx  Sy   Sxx   Syy   Sxy;
    Sx   Sxx Sxy  Sxxx  Sxyy  Sxxy;
    Sy   Sxy Syy  Sxxy  Syyy  Sxyy;
    Sxx Sxxx Sxxy Sxxxx Sxxyy Sxxxy;
    Syy Sxyy Syyy Sxxyy Syyyy Sxyyy;
    Sxy Sxxy Sxyy Sxxxy Sxyyy Sxxyy];

% Construct {k} 
SI = sum(I);
SxI = sum(x.*I);
SyI = sum(y.*I);
SxxI = sum(x.*x.*I);
SyyI = sum(y.*y.*I);
SxyI = sum(x.*y.*I);
k = [SI SxI SyI SxxI SyyI SxyI]';

% Solve
a = C\k;
% Remove background
% First compute background image
[rows, cols] = size(im);
[x, y] = meshgrid( 1:cols, 1:rows );
back = a(1) + a(2)*x + a(3)*y + a(4)*x.*x + a(5)*y.*y +a(6)*x.*y;
im2 = double(im) - back;
im2 = mat2gray(im2); % Convert matrix of type double to image of type double
im2 = im2uint8(im2); % Convert to uint8 image. Although you were not asked
                     % to do this, it is useful to make output same type as
                     % input image.

subplot(1,2,2);
imshow(im2);

end