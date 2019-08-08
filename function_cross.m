function [ Points ] = function_cross( xx,yy,dx,dy,n )
u = linspace(1,n,n)+1;
u = (-1).^u;
u = u*dy;
UX = yy+u;
UY = linspace(xx,xx,n);

Points = [UY;UX];
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


end

