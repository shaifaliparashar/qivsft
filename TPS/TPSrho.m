% TPS kernel
%
% SYNTAX
%   [l,ld]=TPSrho(d)
%  
% INPUT ARGUMENTS
%   -d: a (1xn) array of distances
%
% OUTPUT ARGUMENTS
%   -l: kernelized value
%   -ld: derivative of the kernel
%
% (c) 2013, Daniel Pizarro. dani.pizarro@gmail.com 

% TPS is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
% 
% TPS is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
% for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

function [l,ld]=TPSrho(d)
d2=d;
d2(find(d2==0))=eps;
l=d2.*log(d2);
ld=log(d2)+1;