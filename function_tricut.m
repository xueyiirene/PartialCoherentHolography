function [ V ] = function_tricut( D )
%format A = xa,xb,ya,yb xa<xb ya<yb

LX = D(2)-D(1);
LY = D(4)-D(3);

if LX<LY
    A = [D(1) D(2) D(3)+0.00*(D(4)-D(3)) D(3)+0.36*(D(4)-D(3))];
    B = [D(1) D(2) D(3)+0.30*(D(4)-D(3)) D(3)+0.69*(D(4)-D(3))];
    C = [D(1) D(2) D(3)+0.63*(D(4)-D(3)) D(3)+1.00*(D(4)-D(3))];
    
else
    A = [D(1)+0.00*(D(2)-D(1)) D(1)+0.36*(D(2)-D(1)) D(3) D(4) ];
    B = [D(1)+0.30*(D(2)-D(1)) D(1)+0.69*(D(2)-D(1)) D(3) D(4) ];
    C = [D(1)+0.63*(D(2)-D(1)) D(1)+1.00*(D(2)-D(1)) D(3) D(4) ];
end

V(1,:) = A;
V(2,:) = B;
V(3,:) = C;

end

