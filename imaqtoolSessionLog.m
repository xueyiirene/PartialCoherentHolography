vid = videoinput('pmimaq_2017b', 1, 'PM-Cam 1200x1200');
src = getselectedsource(vid);

vid.FramesPerTrigger = 1;

src.PortSpeedGain = 'Port0-Speed1-100MHz-16bit-Gain1-HDR';

