% 
%   FILE: deliverable_1.m
%   THMMY, 8th semester, Digital Image Process Processing
%   Hough Transform Implementation
%   Author:
%     Moustaklis Apostolos, 9127, amoustakl@auth.gr
%   Hough Transform Implementation
%   Hough Lines Implementation
%   Hough Peaks Implementation

% Load image
I = imread('im2.jpg');
%figure
%imshow(I)

% Convert to grayscale and scale to [0,1]
I = (rgb2gray(I));

global rs;
rs = 0.2;

I=imresize(I,rs);

% Gaussian filter
I = imgaussfilt(I, 4.5);


% Edge filter - use edge()
BW =  edge(I,'Sobel' ); 
figure
imshow(BW)

thetaRes = pi/180;
rhoRes = 1; 
[H,L,res] = myHoughTransform(BW, rhoRes, thetaRes , 12);






function [H,L,res] = myHoughTransform(img_binary, Drho, Drtheta , n )

I= img_binary;
% Rads to degrees
Drtheta = Drtheta * 180/pi;
theBin = -90:Drtheta:90;
% Find the maximum possible d: diagnol length of the image
d = sqrt(size(I,1)^2 + size(I,2)^2);
% Define step size for rhoBin matrix

% Define rho range
rho = -d:Drho:d;
% Rescale rho
rhoBin = 0:Drho:ceil(2*d);
% Initiate hough matrix
H = zeros(length(rhoBin),length(theBin));

%Main Hough Loops 
wb = waitbar(0, 'Computing the Hough Transform');
for i = 1:size(I,1)
    waitbar(i/size(I,1), wb);
    for j = 1:size(I,2)
        if(I(i,j))
            for k = 1:length(theBin) 
                %Calculate R for every theta in the Bin
                tempR = j*cos(theBin(k)* pi/180) + i*sin(theBin(k) * pi/180);
                %Round R to fit the array 
                tempR = round((tempR + d)/Drho)+1;
                H(tempR,k) = H(tempR,k) + 1;
            end
        end
    end
end

close(wb);
 
 imshow(H,[],...
        'XData',theBin,...
        'YData',rho,...
        'InitialMagnification','fit');
 xlabel('\theta (degrees)')
 ylabel('\rho')
 axis on
 axis normal 
 hold on

 
peaks = myHoughPeaks(H,n);

L = peaks ;
%Thresholding the lines next to each other
rhoThres = 0.008*max(rho(:));
thetaThres = 0.03*max(theBin(:));

correctPeaks = 0 ; 
totalIter = 0 ; 


%Thresholding to find the peaks 
    for i = 1 : length(peaks)
         for j = 1 : i
            totalIter = totalIter + 1 ;  
            if (abs(L(j,1) - L(i,1)) < rhoThres )  && ( abs(L(j,1) - L(i,1)) ~= 0 ) 
            L(i,1) = 0;
            elseif (abs(L(j,2) - L(i,2)) < thetaThres )  && ( abs(L(j,1) - L(i,1)) < rhoThres )   && ( abs(L(j,2) - L(i,2)) ~= 0 ) 
            L(i,1) = 0;
            else
               correctPeaks = correctPeaks + 1 ;  
               if correctPeaks == n
                   break;
               end
               continue;
           end
         end       
    end
totalIter  
%Remove from L the thresholded peaks 
L(L(:, 1)== 0, :) = [] ; 
%L(totalIter : end , : ) = [] ;
length(L)
%Remove the Rest 
 L(n+1 : end , : ) = [] ;



 x = theBin(L(:,2));
 y = rho(L(:,1));
 plot(x,y,'s','color','white');


im = imread('im2.jpg');
%Plot the photo and hold on
figure()
imshow(im);
hold on;
global rs;
%Iterate through the peaks array
linePoints = zeros( size(L,1) , 4);
for i = 1:size(L,1)
    %Get the scaled values using theta and rho Scale
    rhoTemp = rho(L(i,1));
    theTemp = theBin(L(i,2));
    if theTemp == 0
        x1 = rhoTemp;
        x2 = rhoTemp;
        y1 = 1;
        y2 = size(I,1);
     
    else
        x1 = 1;
        x2 = size(I,2);
        y1 = (rhoTemp - x1*cos(theTemp*(pi/180))) / sin(theTemp*(pi/180));
        y2 = (rhoTemp - x2*cos(theTemp*(pi/180))) / sin(theTemp*(pi/180));
    end
    %p = [p ; [x1 x2  y1 y2] ];
     plot((1/rs).*[x1,x2],(1/rs).*[y1,y2],'r','LineWidth',2);
     linePoints(i,1) = (1/rs)* x1 ;
     linePoints(i,2) = (1/rs)* x2;
     linePoints(i,3) = (1/rs)* y1;
     linePoints(i,4) = (1/rs)* y2;
   % %scatter([x1,x2],[y1,y2],'*');
    title('Image with hough lines');
end
lineDists = zeros(size(linePoints,1),1);
dist = 0;
for i = 1 : length(linePoints)
    dist = dist + sqrt(  (linePoints(i,1) - linePoints(i,2))^2 +  (linePoints(i,3) - linePoints(i,4))^2 );
end
%Approximately Calculation of the res 
res = size(I,1)*size(I,2) - dist ; 
end




function peaks = myHoughPeaks(H,numPeaks)
    %Double the number of peaks to ensure we detect the correct
    %After the  threshold 
    numPeaks = numPeaks * numPeaks 
    %Copy H into HCopy
    HCopy = H;
    %Initiate peaks matrix
    peaks = zeros(numPeaks,2);
    
    numP = 0;
    %Set a peak threshold 
    threshold = 0.*max(H(:));
    
   while(numP<numPeaks)
       if(threshold<=max(HCopy(:)))
           numP = numP + 1;
           %Find the positions of the max 
           maxP = max(HCopy(:));
           [rows,cols] = find(HCopy == maxP);
           %Make a copy of the cords
           peaks(numP,:) = [rows(1),cols(1)];
           %Set it to zero to find the next
           HCopy(rows(1),cols(1)) = 0;   
       else
           % break if threshold condition is not met
           break;
       end
   end
   % return updated peaks matrix
   peaks = peaks(1:numP,:);
end