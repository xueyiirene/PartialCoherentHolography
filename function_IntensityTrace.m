function [Ica] = function_IntensityTrace(Images,Mask)
% Extract intensity change of a x-y-t image stack in the ROIs of the image. 
% Mask should have the same x-y dimension with Images, but in 2D.
% Mask should be labeled by bwlabel if there are multiple ROIs.
% image is uint12 or uint16, Mask is double. 
% Yi Xue, 2021
N=max(Mask(:));
Mask3D=zeros(size(Mask,1),size(Mask,2),N);
for ii=1:N
    temp=Mask;
    temp(temp~=ii)=0;
    temp(temp~=0)=1;
    Mask3D(:,:,ii)=temp;
end
Ica=zeros(size(Images,3),N);
for ii=1:N
    temp=repmat(Mask3D(:,:,ii),[1,1,size(Images,3)]);
    Ica(:,ii)=squeeze(mean(mean(temp.*double(Images))));
end
end

