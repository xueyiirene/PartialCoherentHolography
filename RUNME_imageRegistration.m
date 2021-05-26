%% image registration
[x,y]=meshgrid(linspace(-400,600,11),linspace(-400,600,11));
DefocusingRadius=zeros(numel(x),1);
TargetRadius=4*abs(x(:).*y(:))/10^5+2;
FinalFrames=function_makespots_ori(Setup,x(:),y(:),DefocusingRadius,TargetRadius);
Mask=FinalFrames(:,:,1);
I=imread('C:\ResearchData\MultiSlice\1-20-20 beadsUSAF\WF1-registration-20xobj\Pos0\img_000000000_Default_000.tif');
%%
registrationEstimator;
cpselect(I,Mask);
%%
t1 = fitgeotrans(movingPoints,fixedPoints,'projective');
Mask_ref = imref2d(size(Mask)); %relate intrinsic and world coordinates
I_registered = imwarp(I,t1,'OutputView',Mask_ref);
figure();
imshowpair(I_registered/max(I_registered(:)),Mask/max(Mask(:)),'blend');
%%
figure();imagesc(I_registered);
threshold=1800;
Mask_temp=I_registered;
Mask_temp(Mask_temp<threshold)=0;
Mask_temp(Mask_temp>=threshold)=1;
Mask_temp=medfilt2(Mask_temp,[5,5]);
FinalFrames=uint8(Mask_temp)*255;

figure();imagesc(FinalFrames);


