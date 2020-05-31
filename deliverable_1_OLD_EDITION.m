% Load image
I = imread('im2.jpg');
figure
imshow(I)

% Convert to grayscale and scale to [0,1]
I = double(rgb2gray(I));
I = I(1:8:end,1:8:end);

%I = medfilt2(I);
%I = imsharpen(I);

%Skip the borders for possible scanning errors
I=I(5:end-5,5:end-5);


% Median blur to smooth the pixture
%I = medfilt2(I);
I = imsharpen(I);

% Edge filter - use edge()
BW = edge(I, 'sobel' );
figure
imshow(BW)

[H,T,R] = myHoughTransform(BW, 1, pi/180 , 5);


function [H,L,res] = myHoughTransform(img_binary, Drho, Drtheta , n )
I = img_binary;
%Im - grayscale image
%Drho - resolution of rhos - scalar
%Drtheta - resolution of theta - scalar



[rows, cols] = size(I);
%Rads to degrees
thetaDeg = Drtheta * 180/pi ; 
theta_maximum = 90;
rho_maximum = (sqrt(rows^2 + cols^2)) ;
thetaScale = -theta_maximum:thetaDeg:theta_maximum - 1;
rhoScale = -rho_maximum:Drho:rho_maximum;

H = zeros(length(rhoScale), length(thetaScale));

wb = waitbar(0, 'Computing the Hough Transform');

for row = 1:rows
    waitbar(row/rows, wb);
    for col = 1:cols
        if I(row, col) > 0
            x = col - 1;
            y = row - 1;
            for theta_ = 1 : length(thetaScale)
                tempR = x*cos(thetaScale(theta_)* pi/180) + y*sin(thetaScale(theta_) * pi/180);
                tempR = round((tempR + rho_maximum)/Drho)+1;
                H(tempR,theta_) = H(tempR,theta_) + 1;
              
            end
        end
    end
end

close(wb);



rhos = zeros(n, 1); % contain rho parameters for lines found in image (indices)
thetas = zeros(n, 1); % contain theta parameters for lines found in image (indices)
H_nms = H; % create copy of accumulator for non-maximal suppression
H_padded= padarray(H_nms, [1, 1], 'replicate'); % take care of boundary problem for NMS by padding with boundary replicates
[rows, cols] = size(H_nms); % size of hough accumulator

for i = 2:rows-1 % to account for padding
    for j = 2:cols-1
        if any(find((H_padded(i-1:i+1, j-1:j+1) > H_padded(i,j)))) > 0 % if any of the neighbors are greater than center pixel
            H_nms(i-1,j-1) = 0; % non-maximal suppression
        end
    end
end

for i = 1:n
   
    
    maxIdx = max(H_nms(:)); % highest score
    [rhoMaxIdx, thetaMaxIdx] = find(H_nms==maxIdx);
    rhos(i) = rhoMaxIdx(1); % add - 1 term to account for padding an extra row and column (3 by 3 filter)
    thetas(i) = thetaMaxIdx(1);
    H_nms(rhoMaxIdx(1), thetaMaxIdx(1)) = 0; % clear the highest scoring cell, then move on
end


L = [rhos , thetas] ; 
res = [];



imshow(H,[],...
       'XData',thetaScale,...
       'YData',rhoScale,...
       'InitialMagnification','fit');
xlabel('\theta (degrees)')
ylabel('\rho')
axis on
axis normal 
hold on

x = thetaScale(L(:,2));
y = rhoScale(L(:,1));
plot(x,y,'s','color','white');



lines = houghlines(I,thetaScale,rhoScale,L,'FillGap',5,'MinLength',7);
length(lines);

figure, imshow(I), hold on

drawHoughLines(I,L , rhoScale , thetaScale);

end



function drawHoughLines(img,peaks,rho,theBin)
%Plot the photo and hold on
figure()
imshow(img);
hold on;
size(peaks,1)
for i = 1:size(peaks,1)
    rhoTemp = rho(peaks(i,1));
    theTemp = theBin(peaks(i,2));
    if theTemp == 0
        x1 = rhoTemp;
        x2 = rhoTemp;
        y1 = 1;
        y2 = size(img,1);
    else
        x1 = 1;
        x2 = size(img,2);
        y1 = (rhoTemp - x1*cos(theTemp*pi/180)) / sin(theTemp*pi/180);
        y2 = (rhoTemp - x2*cos(theTemp*pi/180)) / sin(theTemp*pi/180);
    end
    plot([x1,x2],[y1,y2],'b','LineWidth',2);
    title('Image with hough lines');
    
end

end
