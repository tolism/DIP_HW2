% Load image
I = imread('im3.jpg');
%figure
%imshow(I)

% Convert to grayscale and scale to [0,1]
I = double(rgb2gray(I));
I = I(1:8:end,1:8:end);

%Skip the borders for possible scanning errors
%I=I(5:end-5,5:end-5);

I = medfilt2(I);
%I = imsharpen(I);
%I = imsharpen(I);
% Gaussian filter
%I = imgaussfilt(I, 1);



% Edge filter - use edge()
BW = edge(I, 'Sobel');
%figure
%imshow(BW)

thetaRes = pi/360;
rhoRes = 0.5; 
[H,L,res] = myHoughTransform(BW, rhoRes, thetaRes , 15);
L(:,1) = 8.*L(:,1);
I = imread('im2.jpg');

d = sqrt(size(I,1)^2 + size(I,2)^2);
rho = -d:rhoRes:d;
theBin = -90:thetaRes:90;


%drawHoughLines(I,L,rho,theBin);



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



wb = waitbar(0, 'Computing the Hough Transform');
for i = 1:size(I,1)
    waitbar(i/size(I,1), wb);
    for j = 1:size(I,2)
        if(I(i,j))
            for k = 1:length(theBin) 
                tempR = j*cos(theBin(k)* pi/180) + i*sin(theBin(k) * pi/180);
                tempR = round((tempR + d)/Drho)+1;
                H(tempR,k) = H(tempR,k) + 1;
            end
        end
    end
end

close(wb);
 

% imshow(H,[],...
%        'XData',theBin,...
%        'YData',rho,...
%        'InitialMagnification','fit');
% xlabel('\theta (degrees)')
% ylabel('\rho')
% axis on
% axis normal 
% hold on


peaks = myHoughPeaks(H,n);


L = peaks ;
%Thresholding the lines next to each other
rhoThres = 3;
thetaThres = 2;

 for i = 1 : length(peaks)
     for j = 1 : length(peaks)
          if (abs(L(j,1) - L(i,1)) < rhoThres )  && ( abs(L(j,1) - L(i,1)) ~= 0 ) 
          L(i,1) = 0;
          elseif (abs(L(j,2) - L(i,2)) < thetaThres )  && ( abs(L(j,1) - L(i,1)) < rhoThres )   && ( abs(L(j,2) - L(i,2)) ~= 0 ) 
          L(i,1) = 0;
          else
              continue;
          end
     end
 end

L(L(:, 1)== 0, :) = [] ; 
L ;
res = [];



% x = theBin(L(:,2));
% y = rho(L(:,1));
% plot(x,y,'s','color','white');


drawHoughLines(I,L,rho,theBin);


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


function peaks = myHoughPeaks(H,numPeaks)

% Copy H into HCopy
    HCopy = H;
    %Initiate peaks matrix
    peaks = zeros(numPeaks,2);
    
    numP = 0;
    %Set a peak threshold
    threshold = 0.*max(H(:));
    previousMax = 0;
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