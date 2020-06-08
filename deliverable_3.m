% 
%   FILE: deliverable_1.m
%   THMMY, 8th semester, Digital Image Process Processing
%   Hough Transform Implementation
%   Author:
%     Moustaklis Apostolos, 9127, amoustakl@auth.gr
%   Given the desired angle in rad rotate the input image 


% Load image
I = imread('im2.jpg');
figure
imshow(I)
I=imresize(I,0.2);
%J = imrotate(I,30,'bilinear','loose');
%figure
%imshow(J)


rotImg  =  myImgRotation(I,213*pi/180);
figure
imshow(rotImg);


 
function rotImg = myImgRotation(img , angle)
%Converte the rads to degree
degrees = angle*180/pi;
[m,n,~]=size(img);
%Max diagonal distance 
diagonal = sqrt(m^2+n^2);

wb = waitbar(0, 'Rottating  the image');

for t=1:diagonal
    waitbar(t/diagonal, wb);
   for s=1:diagonal
      i = floor((t-diagonal/2)*cosd(degrees)+(s-diagonal/2)*sind(degrees)+m/2);
      j = floor(-(t-diagonal/2)*sind(degrees)+(s-diagonal/2)*cosd(degrees)+n/2);
      %Check if within image borders
      if i>0 && j>0 && i<=m && j<=n           
         rotImg(t,s,:)=img(i,j,:);
      end
   end
end
close(wb);
 end
  
