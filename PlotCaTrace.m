
I=ToSave.CaTrace;
%%
Rate=50000/80;%NI DAQ frequency 50000Hz, camera frequency 80Hz
t=ToSave.Stim.UUT';
Iref=ToSave.Data;
Iblue=ToSave.Stim.Output(:,1);
Iyellow=ToSave.Stim.Output(:,2);
%%
t_ds=t(1:Rate:end);
%% fluorescence in ROIs
Mask=bwlabel(ToSave.ImageMask(ToSave.Stim.ROIoffset(2):ToSave.Stim.ROIoffset(2)+ToSave.Stim.ROIoffset(4)-1, ...
    ToSave.Stim.ROIoffset(1):ToSave.Stim.ROIoffset(1)+ToSave.Stim.ROIoffset(3)-1),8);
%separate masks into a 3D mask matrix
N=max(Mask(:));
Mask3D=zeros(size(Mask,1),size(Mask,2),N);
for ii=1:N
    temp=Mask;
    temp(temp~=ii)=0;
    temp(temp~=0)=1;
    Mask3D(:,:,ii)=temp;
end
%calcium trace
Ica=zeros(size(I,3),N);
for ii=1:N
    temp=repmat(Mask3D(:,:,ii),[1,1,size(I,3)]);
    Ica(:,ii)=squeeze(mean(mean(temp.*double(I))));
end
%%
figure();plot(t_ds(1:end-15),mean(Ica,2));title('flurescence');
figure();plot(t, medfilt1(Iref,Rate));title('laser photodiode');xlabel('time (s)');
figure();plot(t,Iref);title('raw laser intensity');xlabel('time(s)');
%%
Iref1=medfilt1(Iref,Rate);
Ifluo1=mean(Ica,2);
Ifluo2=Ifluo1(4:end);
Iref1=Iref1(1:Rate:end);
Iref2=Iref1(10:end-9);
%%
figure();plot(t_ds(1:end-18),Ifluo2/max(Ifluo2));title('flurescence measured by camera');
hold on; plot(t_ds(1:end-18),Iref2/max(Iref2));title('laser intensity measured by PD (downsampled)');
xlabel('time (s)');
legend('fluorescence','laser reference');
set(gca,'fontsize',10);
set(gcf,'position',[100,100,400,250]);
axis tight;





%%
Ifluo3=medfilt1(Ifluo2,8);
Iref3=medfilt1(Iref2,8);
%%
figure();plot(Ifluo3./Iref3);
figure();subplot(2, 1, 1);plot(Ifluo3);subplot(2,1,2);plot(Iref3);