function [ Get ] = function_CofC( refAsk,refGet,Ask )
modelterms = [0 0 ; 1 0 ; 0 1 ; 1 1];     %XY spatial calibration model for C_Of_C between SLM and true space
FitX =  polyfitn(refAsk,refGet(:,1),modelterms);
FitY =  polyfitn(refAsk,refGet(:,2),modelterms);

GetX = polyvaln(FitX,Ask);
GetY = polyvaln(FitY,Ask);
Get = [GetX GetY];
end