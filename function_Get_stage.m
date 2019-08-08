function [ positionmicron ] = function_Get_stage( Stage )
position = linspace(0,0,3);
for axisid = 0:2;
toNode1={'0A', '04', hex(fi(int8(axisid))), '00', '00', '00'};
h=hex2dec(toNode1);
fwrite(Stage.serialport,h); 
A = fread(Stage.serialport,12);
sposition = (A(end-3:end));
sposition = typecast(uint8(sposition), 'int32');
position(axisid+1) = sposition;
end

positionmicron = position /Stage.Stepspermicron;
end
