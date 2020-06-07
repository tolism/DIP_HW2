% 
%   FILE: deliverable_1.m
%   THMMY, 8th semester, Digital Image Process Processing
%   Hough Transform Implementation
%   Author:
%     Moustaklis Apostolos, 9127, amoustakl@auth.gr
%   As it's name , just  a scanner for the lazy guys 
%   Change the input picture and the script will find 
%   the number of the photos  
%   Crop them
%   And save them as seperate files 



% Load image
im = imread('im2.jpg');
figure
imshow(im)

% Convert to grayscale and scale to [0,1]
I = (rgb2gray(im));
I=imresize(I,0.2);

% Gaussian filter
I = imgaussfilt(I, 4.6);



% Edge filter - use edge()
BW =  edge(I,'Sobel'); 
figure
imshow(BW)

thetaRes = pi/180;
rhoRes = 1; 
[H,L,res] = myHoughTransform(BW, rhoRes, thetaRes , 12);
%I = imread('im2.jpg');

d = sqrt(size(BW,1)^2 + size(BW,2)^2);

rho = -d:rhoRes:d;
Drtheta = thetaRes * 180/pi;
theBin = -90:Drtheta:90;




%drawHoughLines(I,L,rho,theBin);
%saveas(gcf,'test.png')

%I = imread('test.png');
%I = rgb2gray(I);
corners = myDetectHarrisFeatures( BW);
lineCorn = [];
pointsPerLine = zeros(size(L,1),1) ;
Corners2Line = [];

x = [];
y = [];
linesArray = [];
p = [];

%Calculating the corners that lay at the lines
 for  i = 1:size(L,1)
    
     for j = 1 : size(corners,1)
         rhoTemp = rho(L(i,1));
         theTemp = theBin(L(i,2));
         
         k =  -rhoTemp + corners(j,2)*cos(theTemp*pi/180)  + corners(j,1) *sin(theTemp*pi/180);
     
         Limit = 100 * eps(max(abs([rhoTemp,corners(j,1),corners(j,2),theTemp])));
         if floor(abs(k)) < Limit
             pointsPerLine(i) = pointsPerLine(i) + 1 ;
             c = [corners(j,1) corners(j,2) ];
             lineCorn = [lineCorn; c ];
         end    
     end  
 end

 
 size(lineCorn);
 %Thresholding
 CornerLines = L;
 [rows,cols] = find(pointsPerLine <= 22);
 %To find which lines contain corners
 CornerLines(rows(:),cols(:)) = 0;
 CornerLines(CornerLines(:, 1)== 0, :) = [] ; 

 areParallel =  zeros(size(CornerLines,1) , size(CornerLines,1)) ;
 areVertical =  zeros(size(CornerLines,1) , size(CornerLines,1)) ;
 
 %Find the combinations of the Vertical and Parallel lines
 for i = 1 : length(CornerLines)
     for j = 1 : length(CornerLines)
         %Parallel Check
         if  abs(CornerLines(i,2) - CornerLines(j,2)) < 2  && i ~= j  % && abs(CornerLines(i,1) - CornerLines(j,1)) < 0.48*max(d)
             areParallel(i,j) = 1;
        
         end
          if  abs(CornerLines(i,2) - CornerLines(j,2)) < 92  &&  abs(CornerLines(i,2) - CornerLines(j,2)) > 88  && i ~= j
             areVertical(i,j) = 1;
        
         end
     end
 end
 
 %Isws thn kanw kai auth thn maska me plires diastaseis kai pairnw meta
 parallelCoup = [];
 %Find the Parallel Lines
 for i = 1 : length(CornerLines)
     for j = 1 : length(CornerLines)
         if areParallel(i,j) > 0 && i < j 
             c = [i j];
             parallelCoup = [ parallelCoup ; c ];
         end       
     end
 end
 verticalCoup = zeros(length(parallelCoup) , length(CornerLines)) ;
 
%Find the vertical lines to the parallel's that we calculated
 for i = 1 : length(parallelCoup)
     for j = 1 : length(CornerLines)
         if areVertical(parallelCoup(i,1) , j) > 0  && i ~= j 
            verticalCoup(i,j) = j;
         end    
     end
 end
squareLines = [];
 %Find the pair of parallel / vertical that make a rectangular 
 for i = 1 : length(parallelCoup)
     a = verticalCoup(i,:);
     a = a(a~=0);
     iters = length(a)/2;
     for j = 1 : iters
         a = [ parallelCoup(i,1) parallelCoup(i,2) a(ceil((j+1)/2)) a(ceil((j+2)/2))];
         squareLines = [squareLines ; a];
     end
 end
 
 for i = 1 : size(squareLines,1)
     if squareLines(i,2) == squareLines(i,3)
         squareLines(i,1) =0;
     end
 end
       
%Fixing a bug 
squareLines(squareLines (:, 1)== 0, :) = [] ; 

 %Calculate the Corner Distances 
  harrisCornerDistances = zeros(length(corners),1);
  for i = 1 : length(corners)
       harrisCornerDistances(i) = sqrt( corners(i,1)^2 + corners(i,2));
  end
  
  

 
 %Ean o arithmos den einai 4 tha xreiastei montarisma
 CornersInside = zeros(size(squareLines,1),1);
 for i = 1 : size(squareLines,1)
     %Get the coordinates of the square's corners 
     corn = draw_calculate_interection( squareLines(i,:) , CornerLines , BW ,  rho , theBin );
     if length(corn)~= 0
     if length(corn) == 4 

     % plot(corn(:,1) , corn(:,2) , 'rs' );
     dists = [];
     for k = 1 : length(corn)
       dists =[dists ; sqrt( corn(k,1)^2 + corn(k,2))];
     end
     else
         %Code to process the lines
         display("More than 4 corners ");
         %Will need to do extra process 
     end
     idL = find(dists==min(dists(:)));
     idH = find(dists==max(dists(:)));
      for j = 1 :  size(corners,1) 
          if  corners(j,2) > corn(idL,1) & corners(j,2) < corn(idH,1) & corners(j,1) > corn(idL,2) & corners(j,1) < corn(idH,2) 
                CornersInside(i) = CornersInside(i) + 1 ;
          end
      end
    end
 end

  %To find the non duplicated values   
  stuff = 0 
 [~, ind] = unique(CornersInside(:, 1), 'rows');
 for i = 1 : size(ind)
     if CornersInside(ind(i)) > 0.25*length(corners) 
       stuff = stuff + 1 ;  
     corn = [];
     corn = draw_calculate_interection( squareLines(ind(i),:) , CornerLines , BW ,  rho , theBin   )
     dists = [];
     for k = 1 : length(corn)
       dists =[dists ; sqrt( corn(k,1)^2 + corn(k,2))];
     end
    idL = find(dists==min(dists(:)));
    idH = find(dists==max(dists(:)));
    
    a = im(5*corn(idL,2) : 5*corn(idH,2) , 5*corn(idL,1) : 5*corn(idH,1 ), : );
    figure
    imshow(a);
    hold on
     
     else
         display("Not enough corner information isndie rectangular");
 end
 end
     
 
  

figure
imshow(BW) 
hold on 
plot(lineCorn(:,2) , lineCorn(:,1) , 'rs' );

%Implementation of the function draw_calculate_interactions
%It returns the cordinates of the corners of an image rectangular
 function corn = draw_calculate_interection( squareLines , CornerLines , BW ,  rho , theBin , lineCorn  )
 for i = 1 : size(squareLines,1)
 k = [];
   for j = 1 : size(squareLines,2)
         k =[k ;[CornerLines(squareLines(i,j),1),CornerLines(squareLines(i,j),2)]];   
   end
  %drawHoughLines(BW,k,rho,theBin);
end
      
  
  p = [];
 for i = 1:size(k,1)
    rhoTemp = rho(k(i,1));
    theTemp = theBin(k(i,2));
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
    p = [p ; [x1 x2 y1 y2] ];
    
 end

 lineIntersection = [];
%  length(p)
%        figure
%       title(i);
%       imshow(BW) 
%       hold on 
 for i = 1: length(p)
     for j = 1 : length(p)
         x1 = p(i,1);
         x2 = p(i,2);
         y1 = p(i,3);
         y2 = p(i,4);
         x3 = p(j,1);
         x4 = p(j,2);
         y3 = p(j,3);
         y4 = p(j,4);
         xy = [x1*y2-x2*y1,x3*y4-x4*y3]/[y2-y1,y4-y3;-(x2-x1),-(x4-x3)];
%          scatter(xy(1),xy(2),'*');
         if round(xy(1)) < size(BW,1) && xy(1) > 0 && xy(2) > 0  && round(xy(2)) < size(BW,1) 
  
                     lineIntersection = [lineIntersection ; [ round(xy(1)) round(xy(2))]] ; 
 
         end
              
         end
        
 end
     if length(p) ~= 0
   corn = lineIntersection(1:length(p),:);
      end

 end

 


%Helper function to see if a point is in a specific line
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


%Helper function to find local maxima
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

%Function that detects the Harris Corners
%Returns the cords of the corners in the image I
function corners = myDetectHarrisFeatures(I)

   %The s of the gaussian
   sigma = 1; 
   %Implement the filter
   smoothKernel = fspecial('gaussian',max(1,fix(6*sigma)), sigma);
   %Kerner masks for the gradients 
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
    n = 5000;
   
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
   
  
    corners = [rows,cols]; 
    size(corners)
end 

%Function that performs the Hough Transform
%Returns H array , L with the n peaks and the resolution 
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
rhoThres = 0.006*max(rho(:));
thetaThres = 0.02*max(theBin(:))
%rhoThres = 10;
%thetaThres = 2;

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

%Helper function used to ploot the lines given peaks 
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
    plot([x1,x2],[y1,y2],'r','LineWidth',2);
    %scatter([x1,x2],[y1,y2],'*');
    title('Image with hough lines');
end
hold on;

end

%Helper function to calculate the Hough Transform Peaks
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