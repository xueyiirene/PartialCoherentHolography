function [ Stage ] = function_Start_stage( comport )
try
Stage.Stepspermicron = 40000/8500;    
try
delete(Stage.serialport)    
catch
end
try
clear Stage.serialport    
catch
end
Stage.serialport = serial(comport,'BaudRate',9600);
try
delete(Stage.serialport)    
catch
end
try
clear Stage.serialport    
catch
end
Stage.serialport = serial(comport,'BaudRate',9600);




set(Stage.serialport,'OutputBufferSize',88);


fopen(Stage.serialport);

%disp(['Stage ready on port ' comport]);
catch
    disp('Stage initialization issue');
end
end

