% Load image
I = imread('lena.bmp');
figure
imshow(I)
%I = I(1:8:end,1:8:end);
%J = imrotate(I,30,'bilinear','loose');
%figure
%imshow(J)


rotImg  =  myImgRotation(I,360*pi/180);
figure
imshow(rotImg);


 
function rotImg = myImgRotation(img , angle)
degrees = angle*180/pi;
[m,n,p]=size(img);
mm = m*sqrt(2);
nn = n*sqrt(2);

wb = waitbar(0, 'Rottating  the image');

for t=1:mm
    waitbar(t/mm, wb);
   for s=1:nn
      i = uint16((t-mm/2)*cosd(degrees)+(s-nn/2)*sind(degrees)+m/2);
      j = uint16(-(t-mm/2)*sind(degrees)+(s-nn/2)*cosd(degrees)+n/2);
      if i>0 && j>0 && i<=m && j<=n           
         rotImg(t,s,:)=img(i,j,:);
      end
   end
end
close(wb);
 end
  
