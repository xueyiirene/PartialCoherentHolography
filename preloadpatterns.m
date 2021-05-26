clear all;close all;clc
%%
grid3D=[24,24,2];
N=grid3D(1)*grid3D(2)*grid3D(3);
r=11:-0.1:1;
for j=1:numel(r)
[pts] = round(function_poissonDisc(grid3D,r(j),400));
% A=circshift(pts,[2,1,2]);
% k=1;
% for i=1:52
%     if A(i,1)<=11 && A(i,2)<=11 && A(i,3)<=21
%         B(k,:)=A(i,:);
%         k=k+1;
%     end
% end
if j>1
    pts=cat(1,pts,pts1);
end
temp=sub2ind(grid3D,pts(:,1),pts(:,2),pts(:,3));
[x,y,z] = ind2sub(grid3D,unique(temp,'stable'));
pts1=[x,y,z];
end
%%
temp=sub2ind(grid3D,pts1(:,1),pts1(:,2),pts1(:,3));
C=setdiff(1:N,temp);
C1=C(randperm(length(C)));
[x,y,z]=ind2sub(grid3D,C1');

ptsfinal=cat(1,pts1,[x,y,z]);
%% double check unique
temp=sub2ind(grid3D,ptsfinal(:,1),ptsfinal(:,2),ptsfinal(:,3));
[x,y,z] = ind2sub(grid3D,unique(temp,'stable'));
%%
UXX=linspace(-400,400,grid3D(1));
UZZ=linspace(0,0,grid3D(3));
UX=round(UXX(x));
UY=round(UXX(y));
UZ=round(UZZ(z));
Ind=[x,y,z];
%% 
focicases=[1,4];
seqLcases=[48,48];
%%
figure();
for i=1:121
scatter3(ptsfinal((i-1)*21+1:i*21,1),ptsfinal((i-1)*21+1:i*21,2),ptsfinal((i-1)*21+1:i*21,3),'.');
hold on
end
hold off

%%
for i=1:N
x(i)=Xranges(ptsfinal(i,1));
y(i)=Yranges(ptsfinal(i,2));
z(i)=ZDepths(ptsfinal(i,3));
end
%%
for j=1:121
    DMDFrames_final = gpuArray(false(Setup.DMD.LX,Setup.DMD.LY,210));      
    for i=1:21      
    DMDFrames = function_makespots(Setup,x(i+(j-1)*21),y(i+(j-1)*21),z(i+(j-1)*21),10);   
    DMDFrames_final(:,:,((i-1)*Setup.PointCloud.divider+1):(i*Setup.PointCloud.divider))=DMDFrames;
    end
    save(['Masks11x11x21_10pixelAP\sequence' num2str(j) '.mat'],'DMDFrames_final');
disp(j);
end