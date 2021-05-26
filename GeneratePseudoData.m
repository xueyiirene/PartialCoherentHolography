%generate pseudo data
% Nx=11;
% Ny=11;
% Nz=21;
G=gaussmf(1:numel(Stim.Array),[numel(Stim.Array)/4 numel(Stim.Array)/2]);

B=zeros(size(Stim.Array));
for i = 1:Stim.Npeaks
    B = B+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ).*(1-abs(pp-(Nx+1)/2)/Nx-abs(qq-(Ny+1)/2)/Ny);
end

for pp=1:Nx
    for qq=1:Ny
        A{(pp-1)*Ny+qq,1}=Stim.Output(:,6).*G'.*(1-abs(pp-(Nx+1)/2)/Nx-abs(qq-(Ny+1)/2)/Ny);
    end
end


