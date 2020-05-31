% Load image
%I = imread('im2.jpg');
%figure
%imshow(I)

% Convert to grayscale and scale to [0,1]
%I = rgb2gray(I);
%I = I(1:8:end,1:8:end);
%I = im2double(I / 255);



I  = imread('1.jpg');
I = rgb2gray(I);
%size(I)
corners1 = detectHarrisFeatures(I);
figure
imshow(I);
hold on
plot(corners1);
display('Arithmos corners apo thn ahdia');
size(corners1.Location)
hold off


corners = myDetectHarrisFeatures(I);

figure
imshow(I) 
hold on 
plot(corners(:,1) , corners(:,2) , 'gs' );


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
   threshold = 0.01;
   %Use n if you want to specify the number of points 
   n = 184;
   for i = 2:imgH-1
       for j = 2:imgW-1
           if(isLocalMax(R(i-1:i+1,j-1:j+1)) && R(i,j)> threshold && n > 0  )
               RNonMax(i,j) = 1;
               n = n - 1;
           end
       end
   end
   
   [rows,cols] = find(RNonMax == 1);
   
  
   corners = [cols,rows]; 
    size(corners)

end 


