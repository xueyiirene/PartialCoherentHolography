function [] = function_Stop_stage( Stage )
fclose(Stage.serialport);
%disp('Stage port is closed')
end

