% Load image
I = imread('im2.jpg');
%figure
%imshow(I)

% Convert to grayscale and scale to [0,1]
I = (rgb2gray(I));
I=imresize(I,0.2);

% Gaussian filter
I = imgaussfilt(I, 4.5);



% Edge filter - use edge()
BW =  edge(I,'Sobel'); 
figure
imshow(BW)

thetaRes = pi/180;
rhoRes = 1; 
[H,L,res] = myHoughTransform(BW, rhoRes, thetaRes , 12);
I = imread('im2.jpg');

d = sqrt(size(I,1)^2 + size(I,2)^2);
rho = -d:rhoRes:d;
theBin = -90:thetaRes:90;

pointsArray = [] ;

for i = 1:size(L,1)
    rhoTemp = rho(L(i,1));
    theTemp = theBin(L(i,2));
    if theTemp == 0
        x1 = rhoTemp;
        x2 = rhoTemp;
        y1 = 1;
        y2 = size(BW,1);
    else
        x1 = 1;
        x2 = size(BW,2);
        y1 = (rhoTemp - x1*cos(theTemp*pi/180)) / sin(theTemp*pi/180);
        y2 = (rhoTemp - x2*cos(theTemp*pi/180)) / sin(theTemp*pi/180);
    end
    r = [x1 x2 y1 y2];
    pointsArray = [pointsArray ; r ]; 
end
pointsArray 

%drawHoughLines(I,L,rho,theBin);
saveas(gcf,'test.png')

%I = imread('test.png');
%I = rgb2gray(I);
corners = myDetectHarrisFeatures( BW);
lineCorn = [];

for i = 1 : length(pointsArray)
    for j = 1 : length(corners)
        R = IsPointWithinLine(pointsArray(i,1), pointsArray(i,3), pointsArray(i,2), pointsArray(i,4), corners(j,1), corners(j,2));
        if R
            c = [corners(j,1) corners(j,2) ];
            lineCorn = [lineCorn c ] 
        end
           
            
    end
end
        
% linesEq = zeros(length(pointsArray),1);
% slope = zeros(length(pointsArray) , 1);
% for i = 1 : length(pointsArray)
%     if ( pointsArray(i , 2  ) - pointsArray(i,1 )) ~= 0
%          slope(i) =( pointsArray(i,4 )  - pointsArray(i,3 ) )/( pointsArray(i , 2  ) - pointsArray(i,1 ));
%     
%     else
%         slope(i) = 0;
%         display("Einai katheta to solve ");
%     end
%     for j = 1 : length(corners)
%         g = -corners(j,2) + slope(i)*(corners(j,1) - pointsArray(i,1 ) ) + pointsArray(i,3 ) ;
%         if g == 0 
%             c = [corners(j,1) corners(j,2) ];
%             lineCorn = [lineCorn c ] ;
%         end
%     end
%       
%       
%     
% end
%         
     

figure
imshow(BW) 
hold on 
plot(corners(:,1) , corners(:,2) , 'rs' );

function R = IsPointWithinLine(x1, y1, x2, y2, x3, y3)
% Line equation: y = m*x + b;
Limit = 100 * eps(max(abs([x1,y1,x2,y2,x3,y3])));
if x1 ~= x2
  m   = (y2-y1) / (x2-x1);
  yy3 = m*x3 + y1 - m*x1;
  g = abs(y3 - yy3);
  R   = (abs(y3 - yy3) < 100 * Limit);
else
  R   = (x3 < Limit);
end
end

function flag = isLocalMax(patch)
    pCenter = (size(patch,1)+1)/2;
     [rows,cols] = find(patch == max(patch(:)));
     if(size(rows,1)>1)
         flag = 0;
     elseif(isequal([rows,cols],[pCenter,pCenter]))
         flag = 1;
     else
         flag = 0;
     end
end


function corners = myDetectHarrisFeatures(I)


   sigma = 1; 
   smoothKernel = fspecial('gaussian',max(1,fix(6*sigma)), sigma);
   kerHor = [1 1 1;0 0 0;-1 -1 -1];
   kerVer = [1 0 -1;0 0 0;1 0 -1];
   
   % First derivatives
   Ix = conv2(I,kerVer,'same');
   Iy = conv2(I,kerHor,'same');
   
   % Second degree derivatives
   Ixx = Ix.^2;
   Iyy = Iy.^2;
   Ixy = Ix.*Iy;
   
   % Applying smoothing filter 
   Gxx = conv2(Ixx,smoothKernel,'same');
   Gyy = conv2(Iyy,smoothKernel,'same');
   Gxy = conv2(Ixy,smoothKernel,'same');
   
   % Calculate R matrix
   k = 0.05;
   R = ((Gxx.*Gyy) - (Gxy.^2)) - k * (Gxx+Gyy).^2;
   % Normalize R
   R = R/max(R(:));
   [imgH,imgW] = size(I);
   RNonMax = zeros(imgH,imgW);
   
   %Combinations of thresholding / n gives us the desired corners
   threshold = 0.003;
   %Use n if you want to specify the number of points 
   n = 2532;
   
   wb = waitbar(0, 'Calculating the Corners');
   for i = 2:imgH-1
        waitbar(i/(imgH-1), wb);
       for j = 2:imgW-1
           if(isLocalMax(R(i-1:i+1,j-1:j+1)) && R(i,j)> threshold && n > 0  )
               RNonMax(i,j) = 1;
               n = n - 1;
           end
       end
   end
   close(wb);
   [rows,cols] = find(RNonMax == 1);
   
  
   corners = [cols,rows]; 
    size(corners)
end 





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
rhoThres = 0.008*max(rho(:))
thetaThres = 0.03*max(theBin(:))

correctPeaks = 0 ; 
totalIter = 0 ; 


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
L(L(:, 1)== 0, :) = [] ; 
%L(totalIter : end , : ) = [] ;
length(L)
 L(n+1 : end , : ) = [] ;
res = [];



 x = theBin(L(:,2));
 y = rho(L(:,1));
 plot(x,y,'s','color','white');


drawHoughLines(I,L,rho,theBin);
saveas(gcf,'test.png')

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
    %plot([x1,x2],[y1,y2],'r','LineWidth',2);
    scatter([x1,x2],[y1,y2],'*');
    title('Image with hough lines');
end

end


function peaks = myHoughPeaks(H,numPeaks)

   numPeaks = numPeaks * numPeaks 
% Copy H into HCopy
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

