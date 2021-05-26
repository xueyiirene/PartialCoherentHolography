clear all;close all;clc;

Foldername='SavedData';
Xratio=0.676;%from calibration data
Yratio=0.676;%from calibration data
Zratio=1;
Filename='emx-cre;dio-chrome-s273A_022620a_3_';
cellNum='A';
load([Foldername '/' Filename '.mat']);

% Ind_ori=ToSave.Stim.Ind;
Ind_ori=ToSave.Stim.IndSeq;
Setup.Scorepikes.Method=0;
spikeremoveflag=0;
type=ToSave.type;
StimVoltage=ToSave.Stim.Voltage;

UX=ToSave.Data.UX*Xratio;
UY=ToSave.Data.UY*Yratio;
UZ=ToSave.Data.UZ*Zratio;
Nz=numel(UZ);
NX=numel(UX);
NY=numel(UY);
repeatNum=size(ToSave.Data.XYZ,2);
Stimfreq=ToSave.Stim.FreqHZ;

seqL=ToSave.Stim.Npeaks;
NumofSeq=size(ToSave.Data.XYZ,1);
NumofSpot=ToSave.Stim.NumofSpot;
%%
Score4D=zeros(NX,NY,Nz,repeatNum);
yc=zeros(seqL*NumofSeq,repeatNum);
% switch ToSave.Stim.FreqHZ 
%     case 40
%     DataMatrix3D=zeros(NX*NY,1230,repeatNum);
%     case 50
%     DataMatrix3D=zeros(NX*NY,980,repeatNum);
%     case 80
%     DataMatrix3D=zeros(NumofSeq*seqL,600,repeatNum);
%     case 100
%     DataMatrix3D=zeros(NX*NY,480,repeatNum);
% end
    startpoint=1;
    for j=startpoint:repeatNum
        [cx,cy,cz]=ind2sub([NX,NY,Nz],Ind_ori(:,j));
        for i=1:NumofSeq
                CurrentData=ToSave.Data.XYZ{i,j};
%                 CurrentData=ToSave.Stim.Output(:,1);
                CurrentData=medfilt1(CurrentData,5);
                CurrentData=CurrentData.*ToSave.Stim.CropMask';
                CurrentData(CurrentData==0)=[];
                DataMatrix=reshape(CurrentData,[numel(CurrentData)/seqL,seqL]);
                DataMatrix=DataMatrix';
                DataMatrix3D(1+(i-1)*seqL:i*seqL,:,j)=DataMatrix;
                for p=1:seqL
                    xind=cx(1+(p-1)*NumofSpot+(i-1)*seqL:p*NumofSpot+(i-1)*seqL);
                    yind=cy(1+(p-1)*NumofSpot+(i-1)*seqL:p*NumofSpot+(i-1)*seqL);
                    zind=cz(1+(p-1)*NumofSpot+(i-1)*seqL:p*NumofSpot+(i-1)*seqL);
                        for n=1:NumofSpot
                            [ Score4D(xind(n),yind(n),zind(n),j), odata] = function_score_voltageclamp( Setup,ToSave.Stim,DataMatrix(p,:),spikeremoveflag);  
                        end
                        [yc(p+(i-1)*seqL,j),~]=function_score_voltageclamp( Setup,ToSave.Stim,DataMatrix(p,:),spikeremoveflag);
                end
        end

    end
ScoreVolume=mean(Score4D,numel(size(Score4D)));
figure();
for nz=1:Nz
    imagesc(UX,UY,ScoreVolume(:,:,nz)');
    title([' Z=' num2str(UZ(nz))]);
    caxis([min(ScoreVolume(:)) max(ScoreVolume(:))]);
    axis image;
    xlabel('x(\mum)');ylabel('y(\mum)');
    set(gca,'fontsize',16);   
    saveas(gcf,[Foldername '\' Filename '_plot.fig']);
    saveas(gcf,[Foldername '\' Filename '_plot.tif']);
end
% title(num2str(kk));drawnow;

% figure();set(gcf, 'Position',  [0, 0,1300, 900]);
% for p=1:NX*NY
%     subplot = @(m,n,p) subtightplot (m, n, p);
%     q=sub2ind([NX,NY],ToSave.Ind(p,1),ToSave.Ind(p,2));
%     subplot(NX,NY,q);plot(DataMatrix2D(p,1:20:end));
%     ylim([min(DataMatrix2D(:)) max(DataMatrix2D(:))]);
%     axis off
%     disp(p);
% end
% dim = [0.1 0.66 0.3 0.3];
% str = {['ymin' num2str(min(DataMatrix2D(:))), ' ymax ' num2str(max(DataMatrix2D(:))), '. time(ms) ' num2str(ToSave.Stim.DurationMS)]};
% annotation('textbox',dim,'String',str,'FitBoxToText','on');
% saveas(gcf,[Foldername '\' Filename 'trace' num2str(ProcessData) 'Z' num2str(ToSave.Data.UZ(nz)) '.tif']);
yc=(yc-min(yc(:)))/max(max(yc-min(yc(:))));
figure();imagesc(yc);
xlabel('repeatNum');ylabel('# of illuminations');
set(gcf,'position',[50,50,repeatNum*10,seqL*NumofSeq*8]);  
saveas(gcf,[Foldername '\' Filename '_Yc.fig']);
saveas(gcf,[Foldername '\' Filename '_Yc.tif']);
%%
% C=find(ScoreVolume==repeatNum);
% C1=zeros(NX*NY*Nz,1);
% C1(C)=1;

h=zeros(NX*NY*Nz/NumofSpot*repeatNum,NX*NY*Nz);
Ind=Ind_ori(:);
for i=1:size(h,1)   
    tempInd=Ind((1+(i-1)*NumofSpot):i*NumofSpot);
    h(i,tempInd)=1;
end
figure();imagesc(h);
xlabel('all pixels');ylabel('# of illuminations');
set(gcf,'position',[50,50,NX*NY*Nz,seqL*NumofSeq*8]);  
saveas(gcf,[Foldername '\' Filename '_h.fig']);
saveas(gcf,[Foldername '\' Filename '_h.tif']);
% S=20*h\ScoreVolume(:);
% nnz(S)

%% generate simulation data to test optimization algorithm
x_gt=zeros(NX,NY,Nz);
x_gt(11:13,11:13,5:7)=4;
x_gt(10:14,10:14,6)=4;

for k=1:repeatNum
yc(:,k)=h*x_gt(:)+randn(NumofSeq*seqL,1);
% yc(:,k)=h*x_gt(:);
end
figure();imagesc(yc);
xlabel('repeatNum');ylabel('# of illuminations');
%% test how different number of repetition affects results
% yc_ori=yc;h_ori=h;
yc=yc_ori(:,1:5);
h=h_ori(1:80*5,:);

%% Optimization
% h=h';
% x0=rand(NX*NY*Nz,1);
tic
reg='TV';
rng('shuffle');
x0=rand(NX*NY*Nz,1);
h_T=conj(h');
% Lmax=abs(max(eig(A_T*HStack)));% 0<L<Lmax
L=1500;
alpha=0.01;
lam=2;
max_iter = 500;
xk=x0;
tk=1;%initial tk
u=xk;
Yc=yc(:);
    for i=1:max_iter
        Uk=u;
        residue=h*Uk-Yc;
        gradf = h_T*residue;
        Uk1 = Uk - 1/L*gradf;
        switch reg 
            case 'sparse'
                xk1= wthresh(Uk1,'s',lam/L);
            case 'TV'
                Uk13D=reshape(Uk1,[NX,NY,Nz]);
                K=3*Uk13D-circshift(Uk13D,1,1)-circshift(Uk13D,1,2)-circshift(Uk13D,1,3);
                xk1=wthresh(Uk1-alpha*K(:),'s',lam/L);
        end
        tk1=1+sqrt(1+4*tk^2)/2;
        u=xk1+(tk-1)/tk1*(xk1-xk);
        xk=xk1;
%         figure(2);imagesc(xk);
        if i==1
            loss=sqrt(sum((h*xk-Yc).^2))+lam*sum(abs(xk(:)));
        else
            loss=cat(1,loss,sqrt(sum((h*xk-Yc).^2))+lam*sum(abs(xk(:))));
%             loss=cat(1,loss,sqrt(sum(sum((xk-S.*x).^2))));
        end
        if i>2 && loss(end)>loss(end-1)
            disp('cost increase');
            break;
%         elseif sum(abs(xk(:)-Sx(:)))<0.001
%             disp('error is too small');
%             break;
        else
%             figure(1);plot(loss);drawnow;
        end
        
        if rem(i,10)==0
            x_obj=reshape(xk,[NX,NY,Nz]);
            figure(3);
            for nz=1:Nz
                imagesc(UX,UY,x_obj(:,:,nz)');drawnow;
                title([' i=' num2str(i)]);
                caxis([min(x_obj(:)) max(x_obj(:))]);
%                 axis image;
            end
        end
    end
xfinal=reshape(xk,[NX,NY,Nz]);   
figure();
for nz=1:Nz
imagesc(UX,UY,xfinal(:,:,nz)');
title([' Z=' num2str(UZ(nz))]);
caxis([min(xfinal(:)) max(xfinal(:))]);
axis image;
xlabel('x(\mum)');ylabel('y(\mum)');
set(gca,'fontsize',16);   
saveas(gcf,[Foldername '\' Filename '_x_recon.fig']);
saveas(gcf,[Foldername '\' Filename '_x_recon.tif']);
end

toc
figure();plot(loss);title(['FISTA+' reg]);xlabel('iteration #');ylabel('loss');
dim = [0.5 0.55 0.3 0.3];
switch reg
    case 'TV'
        str = {['L= ' num2str(L), ', \alpha= ' num2str(alpha), ', lam= ' num2str(lam)]};
    case 'sparse'
        str = {['L= ' num2str(L), ', lam= ' num2str(lam)]};
end
annotation('textbox',dim,'String',str,'FitBoxToText','on');
saveas(gcf,[Foldername '\' Filename '_loss.fig']);
saveas(gcf,[Foldername '\' Filename '_loss.tif']);

%% hot plot
% figure();
% set(gcf,'position',[100, 100, 500, 400]);
% % subplot(1,2,1);
% % imagesc(UX,UY,Score3D_random(:,:,5)');
% % axis image;xlabel('x(\mum)');ylabel('y(\mum)');
% % caxis([min(min(Score3D(:)),min(Score3D_random(:))),max(max(Score3D(:)),max(Score3D_random(:)))]);
% % colormap('hot');
% % set(gca,'fontsize',16);
% % subplot(1,2,2);
% % imagesc(UX,UY,Score3D(:,:,5)');
% % axis image;
% % caxis([min(min(Score3D(:)),min(Score3D_random(:))),max(max(Score3D(:)),max(Score3D_random(:)))]);
% % colormap('hot');xlabel('x(\mum)');ylabel('y(\mum)');
% % set(gca,'fontsize',16);
% 
% for nz=1:Nz
% subplot(1,Nz,nz);
% imagesc(UX,UY,S(:,:,nz)');
% title(['Z = ' num2str(ToSave.Data.UZ(nz),2) '\mum']);
% axis image;
% xlabel('x(\mum)');ylabel('y(\mum)');
% set(gca,'fontsize',16);
% axis off;
% caxis([min((S(:))) max((S(:)))]);colormap('hot');
% end
% saveas(gcf,[Foldername '\' Filename '_hotmap.fig']);
% saveas(gcf,[Foldername '\' Filename '_hotmap.pdf']);
%%
figure(6);
for nz=1:Nz
    subplot(2,4,nz);
    imagesc(UX,UY,ScoreVolume(:,:,nz)');
    title([' Z=' num2str(UZ(nz))]);
    caxis([min(ScoreVolume(:)) max(ScoreVolume(:))]);
    axis image;
end

% xlabel('X (\mum)');ylabel('Y (\mum)');
% set(gca,'FontSize',16);
%%
saveas(gcf,[Foldername '/' Filename 'Volume.tif']);
saveas(gcf,[Foldername '/' Filename 'Volume.fig']);
save([Foldername '/Score4D_' Filename '_' cellNum '.mat'],'Score4D');
%%
artifact=zeros(24,24,8);
for i=1:24
     cx=Ind(1+(i-1)*seqL:NumofSpot+(i-1)*seqL,1);
     cy=Ind(1+(i-1)*seqL:NumofSpot+(i-1)*seqL,2);
     cz=Ind(1+(i-1)*seqL:NumofSpot+(i-1)*seqL,3);
    artifact(cx,cy,cz)=1;
end
figure(3);
for nz=1:Nz
    subplot(2,4,nz);
    imagesc(UX,UY,artifact(:,:,nz)');
    title([' Z=' num2str(ToSave.Data.UZ(nz),2)]);
    caxis([min(artifact(:)) max(artifact(:))]);
    axis image;
end
ScoreVolume=ScoreVolume-artifact*repeatNum;
%% temporal sequence


%% plot part of the region
Indrandom=sub2ind([NX,NY],ToSave.Ind(:,1),ToSave.Ind(:,2));
[indx_roi,indy_roi]=meshgrid(1:20,1:20);
Indaim=sub2ind([NX,NY],indx_roi(:),indy_roi(:));
[~,irandom,iroi] = intersect(Indrandom,Indaim,'sorted'); 
figure();set(gcf, 'Position',  [100, 100,size(indx_roi,2)*50, size(indx_roi,1)*50]);
for p=1:numel(Indaim)
    subplot = @(m,n,p) subtightplot (m, n, p);
subplot(size(indx_roi,1),size(indx_roi,2),p);
plot(DataMatrix2D(irandom(p),1:40:end));
% ylim([-0.8 -0.07]);
ylim([min(min(DataMatrix2D(irandom,:))) max(max(DataMatrix2D(irandom,:)))]);
axis off
disp(p);
end
dim = [0.1 0.66 0.3 0.3];
str = {['ymin' num2str(min(min(DataMatrix2D(irandom,:)))), ' ymax ' num2str(max(max(DataMatrix2D(irandom,:)))), '. time(ms) ' num2str(ToSave.Stim.DurationMS)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
saveas(gcf,[Foldername '\' Filename 'trace_roi1_Z' num2str(ToSave.Data.UZ(nz)) '.fig']);
disp([min(min(DataMatrix2D(irandom,:))) max(max(DataMatrix2D(irandom,:)))]);
%% plot seal test to check the condition of cell
sealtest=zeros(size(ToSave.Data.XY));
for j=1:repeatNum
    for nz=1:Nz
        for i=1:NumofSeq
            if isempty(ToSave.Data.XY{i,j,nz})
                continue;
            else
                temp=ToSave.Data.XY{i,j,nz};
                sealtest(i,j,nz)=max(temp(1:2500));
            end
        end
    end
end
sealtest1=permute(sealtest,[1,3,2]);
figure();plot(sealtest1(:));
axis tight;grid on; 
ylabel('mV');
xlabel('Seal test of each sequence');
saveas(gcf,[Foldername '\' Filename 'sealtest.fig']);


%%
% if ProcessData==1
%     SXY2D=reshape(mean(SXY,2),[numel(UX) numel(UY)]);
% else
%     SXY2D=mean(SXY,3);
% end
%%
% [~,Ind]=max(mean(SXY,2));
% if numel(Ind)>1
%     [~,ind]=min(abs(Ind-length(SXY)/2));
%     Ind=Ind(ind);
% end
% 
% [I,J]=ind2sub([length(UX) length(UY)],Ind);
%%
%  figure();imagesc(log(Score2D'));
row=29;col=20;startpoint=2;
Indrandom=sub2ind([NX,NY],ToSave.Ind(:,1),ToSave.Ind(:,2));
q=sub2ind([NX,NY],row,col);
[v]=find(Indrandom==q);
A=DataMatrix3D(v,:,:);
A=squeeze(A);
if startpoint==1
    A=A';
end
t=0:ToSave.Stim.UT(2):ToSave.Stim.UT(2)*(size(A,1)-1);
[~,peak]=max(abs(mean(A(:,startpoint:end),2)));
a1=mean(A(peak:end,startpoint:end),2);
t1=t(peak:end)';
[fitresult, gof] = function_createMonoExpFit(t1, abs(a1));
temp=coeffvalues(fitresult);tau=-1000/temp(2);%exponential decay time

figure();
for i=startpoint:repeatNum
    plot(t,A(:,i),'color',[1-i/repeatNum, 0, i/repeatNum]);hold on;
end
hold off
dim = [0.7 0.6 0.3 0.3];
str = {['\tau=' num2str(tau) 'ms']};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
xlim([min(t) max(t)]);
xlabel('time (s)');ylabel('photocurrent (nA)');
grid on;
set(gca,'FontSize',16);
saveas(gcf,[Foldername '\' Filename 'one_trace_' num2str(row) '_' num2str(col) '_Z' num2str(ToSave.Data.UZ(nz)) 'V2.6.tif']);
 %% plot decay time map
 DecayTimeMap=zeros(NX,NY);
 startpoint=2;
 
 for yy=1:NY
     for xx=1:NX
        row=xx;col=yy;
        Indrandom=sub2ind([NX,NY],ToSave.Ind(:,1),ToSave.Ind(:,2));
        q=sub2ind([NX,NY],row,col);
        [v]=find(Indrandom==q);
        A=DataMatrix3D(v,:,:);
        A=squeeze(A);
        if startpoint==1
            A=A';
        end
        [peakVal,peak]=max(abs(mean(A(:,startpoint:end),2)));
        baseline=A(1:10,startpoint:end);
        if peakVal<abs(mean(baseline(:)))+0.08
            tau=nan;
        else
            a1=mean(A(peak:end,startpoint:end),2);
            t=0:ToSave.Stim.UT(2):ToSave.Stim.UT(2)*(size(A,1)-1);
            t1=t(peak:end)';
            [fitresult, gof] = function_createMonoExpFit(t1, abs(a1));
            temp=coeffvalues(fitresult);tau=-1000/temp(2);%exponential decay time
        end
        DecayTimeMap(col,row)=tau;
     end
 end
 figure();
 imagesc(UX,UY,DecayTimeMap);
axis image;xlabel('X (\mum)');ylabel('Y (\mum)');
colorbar;
set(gca,'FontSize',16);
title(['Decay time map (ms), Z=' num2str(ToSave.Data.UZ(nz))]);
saveas(gcf,[Foldername '/' Filename 'DecayTimePlot_Z' num2str(ToSave.Data.UZ(nz)) '.tif']);
saveas(gcf,[Foldername '/' Filename 'DecayTimePlot_Z' num2str(ToSave.Data.UZ(nz)) '.fig']);
