function Setup = function_initializeDAQ(varargin)
Setup.Daq = daq.createSession('ni');
Setup.Daq.Rate = 50000;
addAnalogOutputChannel(Setup.Daq,'Dev2','ao0','Voltage');
addAnalogOutputChannel(Setup.Daq,'Dev2','ao1','Voltage');
addAnalogOutputChannel(Setup.Daq,'Dev2','ao2','Voltage');
addAnalogOutputChannel(Setup.Daq,'Dev2','ao3','Voltage');
addDigitalChannel(Setup.Daq,'Dev2','Port0/Line0:1','OutputOnly');
end

