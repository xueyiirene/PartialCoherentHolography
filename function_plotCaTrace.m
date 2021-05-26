function [CaTrace, dF, Mask] = function_plotCaTrace(Images, Stim)
%plot dF/F trace from the input x-y-t images
% Yi Xue, 3/26/2021
I1=double(Images); %kick out bad images
background=repmat(mean(mean(I1(1:50,1:50,:))),[size(I1,1), size(I1,2),1]);
I1=I1-background;
% figure();imshow3D(I1);
% generate Mask
Imean=uint16(mean(I1,3));

mask1=imbinarize(Imean,'global');
mask2=imbinarize(Imean,'adaptive', 'Sensitivity', 0.57);
Mask=medfilt2(mask1.*mask2,[15, 15]);
Mask(1:70,:)=0;
Mask(140:end,:)=0;
Mask(:,1:70)=0;
Mask(:,140:end)=0;

figure(256);imagesc(Imean);axis image;colormap('gray');caxis([min(Imean(:)) max(Imean(:))]);colorbar;
green = cat(3, zeros(size(Imean)),ones(size(Imean)),zeros(size(Imean)));
hold on
h=imshow(green);
hold off
set(h,'AlphaData',Mask*0.1);
drawnow;

% Intensity trace
Itemp=I1.*Mask;
CaTrace=squeeze(sum(sum(Itemp)))/nnz(Mask);
% dF/F trace
dF=zeros(size(CaTrace));
w=round((1000/Stim.Exposure)/Stim.cutoffFreq*0.4);
for ii=1:numel(CaTrace)
    if ii<=w
        F0=min(CaTrace(1:ii+w));
    elseif ii>numel(CaTrace)-w
        F0=min(CaTrace(ii-w:end));
    else
        F0=min(CaTrace(ii-w:ii+w));
    end
    dF(ii)=(CaTrace(ii)-F0)/F0;
end
end


