% TPS Warp Interpolation of points, giving derivatives
%
% SYNTAX
%   [C]=TPSGenerateCenters(nC,KLims,rnd)
%  
% INPUT ARGUMENTS
%   -nC: the warp will use nC^r centers
%   -KLims: Bounds to spread the centers. For r=1 KLims=[umin,umax]. For r=2
%   KLims=[umin,umax,vmin,vmax]. For r=3
%   KLims=[umin,umax,vmin,vmax,zmin,zmax]. r is the dimmension of input
%   points
%   -rnd =0 uses a regular grid for the control centers
%   -rnd =1 spreads the control centers randomly   
%
% OUTPUT ARGUMENTS
%   %   -C: TPS centers in the template
%
% (c) 2013, Daniel Pizarro. dani.pizarro@gmail.com 
%
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


function [C]=TPSGenerateCenters(nC,KLims,rnd)
if(nargin<3)
    rnd=0; % regular grid
end

switch(length(KLims))
    case 2
        if(rnd==0)
        C=linspace(KLims(1),KLims(2),nC)';
        else
        C=(KLims(2)-KLims(1))*rand(nC,1)+KLims(1);    
        end
        
    case 4
        if(rnd==0)
        [Cx,Cy]=meshgrid(linspace(KLims(1),KLims(2),nC),linspace(KLims(3),KLims(4),nC));
        C=[Cx(:),Cy(:)];
        else
        Cx=(KLims(2)-KLims(1))*rand(nC,nC)+KLims(1);
        Cy=(KLims(4)-KLims(3))*rand(nC,nC)+KLims(3);
        C=[Cx(:),Cy(:)];
        end       
    case 6
        if(rnd==0)
             [Cx,Cy,Cz]=meshgrid(linspace(KLims(1),KLims(2),nC),linspace(KLims(3),KLims(4),nC),linspace(KLims(5),KLims(6),nC));
        C=[Cx(:),Cy(:),Cz(:)];
        else
            Cx=(KLims(2)-KLims(1))*rand(nC,nC,nC)+KLims(1);
        Cy=(KLims(4)-KLims(3))*rand(nC,nC,nC)+KLims(3);
        Cz=(KLims(6)-KLims(5))*rand(nC,nC,nC)+KLims(5);
        C=[Cx(:),Cy(:),Cz(:)];
        end
        
    otherwise
    C=[];
    disp('Error: KLims must be a 1x2 [umin,umax] 1x4 [umin,umax,vmin,vmax] or 1x6 [xmin,xmax,ymin,ymax,zmin,zmax] vector')
end

end


