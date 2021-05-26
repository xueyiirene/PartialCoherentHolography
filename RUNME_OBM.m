%%
foldername1='WF9';
foldername2='WF10';
Iraw1=imread(['C:\ResearchData\PCH V2\ControlCode\PartialCoherentHolography-master\IR OBM\' foldername1 '\Pos0\img_000000000_Default_000.tif']);
Iraw2=imread(['C:\ResearchData\PCH V2\ControlCode\PartialCoherentHolography-master\IR OBM\' foldername2 '\Pos0\img_000000000_Default_000.tif']);
% BG=imread('C:\ResearchData\PCH V2\ControlCode\PartialCoherentHolography-master\IR OBM\IRBG.tif');

% I1=double(Iraw1)-double(BG);
% I2=double(Iraw2)-double(BG);
I1=double(Iraw1);
I2=double(Iraw2);
% I1=I1./medfilt2(I1,[300 300]);
% I2=I1./medfilt2(I2,[300 300]);

I1=I1./imgaussfilt(I1,100);
I2=I2./imgaussfilt(I2,100);
I1=(I1-min(I1(:)))/max(max(I1-min(I1(:))));
I2=(I2-min(I2(:)))/max(max(I2-min(I2(:))));
Ifinal=I1-I2;
figure();subplot(1,3,1);imagesc(I1);axis image;colormap('gray');
subplot(1,3,2);imagesc(I2);axis image;colormap('gray');
subplot(1,3,3);imagesc(Ifinal);axis image;colormap('gray');caxis([-0.1 0.2]);
Ig=edge(I2);
figure();imagesc(Ig);
