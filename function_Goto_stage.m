function [ ] = function_Goto_stage( Stage,positionmicron )
position = round(positionmicron * Stage.Stepspermicron);


for axisid = 0:2;
mot = hex(fi(int32(position(axisid+1))));
toNode1={'53', '04', '06', '00', '00', '00', hex(fi(int8(axisid))), '00', mot(7:8), mot(5:6), mot(3:4), mot(1:2)};
h=hex2dec(toNode1);
fwrite(Stage.serialport,h); 
end

end

