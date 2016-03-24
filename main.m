% add path to where the data are stored
addpath('/Users/paul5/Google Drive/arabidopsis/automated quantification test')
% importing images
image = imread('/Users/paul5/Google Drive/arabidopsis/automated quantification test/1 COL/Series008_z0.tif');
image = image(1:1000,:,1);

%%

figure, subplot(1,2,1), imagesc(image), subplot(1,2,2), hist(double(image(:)))

%%
imagetemp = zeros(size(image));
for i = 1:size(image,1),
    for j = 1:size(image,2),
        if image(i,j)>25,
            imagetemp(i,j) = image(i,j);
        end
    end
end

figure, imagesc(imagetemp)

%% edge detection

[~, threshold] = edge(image, 'sobel');
fudgeFactor = .7;
BWs = edge(image,'sobel', threshold * fudgeFactor);
figure, imshow(BWs), title('binary gradient mask');

se90 = strel('line', 3,90);
se0 = strel('line', 3, 0);

BWsdil = imdilate(BWs, strel(ones(6)));%[se90 se0]);
figure,subplot(1,2,1), imshow(BWsdil), title('dilated gradient mask');
subplot(1,2,2),imshow(image),

%% segmentation

L = bwlabel(1-BWsdil);

figure, hist(L(:));
figure, imagesc(L);

%%

L_filtered = L;
I = 1:max(max(L));
c = 1;
for i = I;
    indtemp = (find(L(:) == i)) ;
    if length(indtemp) < 250,
        L_filtered(indtemp) = 0;
    else
        L_filtered(indtemp) = c;
        c = c+1;
    end
end
figure, imagesc(L_filtered);

%% looking at all cells

L_filtered_2 = zeros(size(L_filtered));

figure,
for i = 1:max(max(L_filtered));
    mask = (find(L_filtered(:) == i));
    L_filtered_2(mask) = max(max(L_filtered))-i;
    L_filtered_2(~mask) = 0;
    imagesc(L_filtered_2);
    title(i);
    pause(0.25);
end


