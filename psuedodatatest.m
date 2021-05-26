x0=zeros(20,20);
x0(10:14,4:10)=1;
y=zeros(100*Stim.repeatNum,1);
h=zeros(20,20,100*Stim.repeatNum);
%%

for j=1:Stim.repeatNum
    Stim.IndSeq(:,j)=randperm(Nx*Ny);
    for i=1:100
        tempInd=Stim.IndSeq((1+(k-1)*seqL+(i-1)*Stim.NumofSpot):((k-1)*seqL+i*Stim.NumofSpot),j);
        for n=1:Stim.NumofSpot
            h(xind(tempInd(n)),yind(tempInd(n)),i+(j-1)*100)=1;
        end
        y(i+(j-1)*100)=sum(sum(h(:,:,i+(j-1)*100).*x0));
    end
end
h1=reshape(h,[1000,400]);
x=reshape(h1\y,20,20);

i=1;

for j=1:10
Score3D=zeros(Nx,Ny,Nz);
for k=1:seqL
            tempInd=Stim.IndSeq((1+(i-1)*seqL+(k-1)*Stim.NumofSpot):((i-1)*seqL+k*Stim.NumofSpot),j);
            for n=1:Stim.NumofSpot
                Score3D(xind(tempInd(n)),yind(tempInd(n)),nz) = y(k+(j-1)*100);
            end
end
Score(:,:,j)=Score3D;
figure();subplot(1,2,1);imagesc(Score3D);subplot(1,2,2);imagesc(sum(Score,3));drawnow;
end
%%
x0=zeros(1600,2560);
x0(560:900,1000:1390)=1;
%%
for zz=1:seqL
    y(zz)=sum(sum(x0.*DMDFrames_final(:,:,zz)));
end


for i=1:NumofSeq
        
        for k=1:seqL
            tempInd=Stim.IndSeq((1+(i-1)*seqL+(k-1)*Stim.NumofSpot):((i-1)*seqL+k*Stim.NumofSpot),j);
            for n=1:Stim.NumofSpot
              Score3D(xind(tempInd(n)),yind(tempInd(n)),nz)= y(k);
            end
        end
end
Score(:,:,j)=Score3D;
figure();subplot(1,2,1);imagesc(Score3D);subplot(1,2,2);imagesc(sum(Score,3));

