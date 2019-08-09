function [right]=binarysearch(nx,x,value,d)

% Given an array and a value, returns the index of the element that
% is closest to, but less than, the given value.
% Uses a binary search algorithm.
% "delta" is the tolerance used to determine if two values are equal
% if ( abs(x1 - x2) <= delta) then
% assume x1 = x2
% endif

% example d = 1e-9; nx=length(atmo.x); x=atmo.x; value g.x(i)

left = 1;
right = nx;

while (left < right)
%    if (left > right) then
%            exit
%       endif
  middle = round((left+right) / 2.0);
  if ( abs(x(middle) - value) <= d)
      binarySearch = middle;
      return
 elseif (x(middle) > value) 
     right = middle - 1;
 else
     left = middle + 1;
 end
end
        
binarysearch = right;


