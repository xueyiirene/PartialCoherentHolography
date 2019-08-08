function [  ] = function_livedisplay( f,Setup,Result,i ,j)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
subplot(2,2,1)
plot(Setup.UT,Result.Signal{i}); hold on
a = Result.peakdip{i};
scatter(a(1), a(2),'blue','filled')
scatter(Setup.UT(Setup.UT<Setup.Ephys.PatchTest.BeginPause),Result.Baseline{i},'red','filled')
xlabel('Time(s)'); ylabel('Patch Signal (Volts)');
hold off
select = Result.Done==1;
subplot(2,2,2)
scatter3(Result.StagePositions(3,:),Result.StagePositions(2,:),Result.StagePositions(1,:), 'blue')
xlabel('X \mum');ylabel('Y \mum');zlabel('Z \mum'); hold on;
scatter3(Result.StagePositions(3,i),Result.StagePositions(2,i),Result.StagePositions(1,i),'filled', 'red');hold on;
axis([min(Result.StagePositions(3,:)-1) max(Result.StagePositions(3,:)+1) min(Result.StagePositions(2,:)-1) max(Result.StagePositions(2,:)+1) min(Result.StagePositions(1,:)-1) max(Result.StagePositions(1,:)+1)  ])
scatter3(Result.StagePositions(3,select),Result.StagePositions(2,select),Result.StagePositions(1,select),[],Result.PPSF(select,j),'filled')
xlabel('X \mum');ylabel('Y \mum');zlabel('Z \mum'); colorbar; title('Peak Current (measured Volts)'); hold off;
view(3)
subplot(2,2,3)
plot(Setup.UT,Result.LightOutputs(j,:));
xlabel('Time(s)'); ylabel('Light Signal (Volts)')

drawnow;

end

