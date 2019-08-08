function [ ] = function_Zero_stage( Stage )
for axisid = 0:2;
toNode1={'09', '04', '06', '00', '00', '00', hex(fi(int8(axisid))), '00', '00', '00', '00', '00'};
h=hex2dec(toNode1);
fwrite(Stage.serialport,h); 
end
end

