function [ logvec ] = function_logvec( range,n )
if n == 0; logvec = zeros(0);
elseif n==1; logvec = 0;
else    
vec = linspace(-log(range),log(range),n);
signvec = sign(vec);
vec = exp(abs(vec));
logvec = vec.*signvec;
end
end

