function [num] = count_binary_beads_for_pipe(im_path)

%% binarizing bead images



I = tiffreadVolume(im_path); % laod z stack

I = max(I, [], 3); % get max projection

% figure
% imagesc(I);

I2 = imgaussfilt(I, 40);

% I2 = medfilt2(I, [100,100]);

% I = imadjust(I);

% I2 = reshape(I, [] , 1);
% 
% lv = median(I);
% 
% I = I - lv;

I = imsubtract(I, I2);

% figure 
% imagesc(I);


%%

Igray = mat2gray(I);

% imshow(Igray);


%%
%convert into grayscale image. 
k = Igray;
%calculate threshold using Otsu's method. 
level= graythresh(k); 
  
%convert into binary image using level. 
k1=imbinarize(k,level); 
  
%display the binarized image. 
% imagesc(k1); 


[L, num] = bwlabel(k1);

fprintf('Found approximately %d beads', num)

end
