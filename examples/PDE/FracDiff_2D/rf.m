function [rhs] = rf(f,t,J,s)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
rhs=[];
for ns=1:s
    rhs=[rhs reshape(f(t(ns)),(J-1)^2,1)];
end
end

