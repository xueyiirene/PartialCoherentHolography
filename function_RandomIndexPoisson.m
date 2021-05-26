function [ptsfinal] = function_RandomIndexPoisson(Nx,Ny)
%for a 1:Nx*Ny grid point, generate a random sequence of it using poisson
%disc
spacing=4;%default spacing
nPts0=round(Nx*Ny/5);
nPts=nPts0;
[pts] = round(function_poissonDisc([Nx Ny],spacing,nPts,0));
pts1=pts;
for spacing=3.5:-0.5:1.5
    nPts=round(max((Nx*Ny-numel(pts1))/2,nPts0));
    [pts] = round(function_poissonDisc([Nx Ny],spacing,nPts,0));
    if isempty(pts)
        continue;
    else
        pts=cat(1,pts,pts1);
        temp=sub2ind([Nx,Ny],pts(:,1),pts(:,2));
        [x,y] = ind2sub([Nx,Ny],unique(temp,'stable'));
        pts1=[x,y];
    end
end
temp=sub2ind([Nx,Ny],pts1(:,1),pts1(:,2));
C=setdiff(1:Nx*Ny,temp);
C1=C(randperm(length(C)));
[x,y]=ind2sub([Nx,Ny],C1');
pts2=cat(1,pts1,[x,y]);

Dist=sqrt(diff(pts2(:,1)).^2+diff(pts2(:,2)).^2);
Ind=find(Dist<3);
Ar=pts2(Ind,1);
Ac=pts2(Ind,2);

for i=1:round(numel(Ind)/2)
Ar([i numel(Ind)-i+1])=Ar([numel(Ind)-i+1 i]);
Ac([i numel(Ind)-i+1])=Ac([numel(Ind)-i+1 i]);
end
pts2(Ind,1)=Ar;
pts2(Ind,2)=Ac;
ptsfinal=pts2;
end

