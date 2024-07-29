function [sag] = calc_sag(R, N_pixels_obj, N_pixels_mla, interval, lenscentre)
% This function calculates the sagitta (sag) value of a lens. 

% INPUTS
% R - radius of curvature from the lens maker's equation (m)
% interval - can be the manufacturing resolution or camera pixel size (m)
% lenscentre - the xy coordinates of MLA centres (pixels)

d = ( -N_pixels_obj / 2 ):( N_pixels_obj / 2 - 1 );
d = d * interval;
lenscentre = lenscentre * interval;

s = R - sqrt(R^2 - ( ( ( 2 * N_pixels_mla / sqrt(3) ) ) / 2 * interval ) .^ 2 ); %take diagonal to make sag extend over whole hexagon area

% tilt in X direction
[x,y] = meshgrid(d,d);
rho = sqrt( ( ( x-lenscentre(1) ) .^ 2 ) + ( ( y - lenscentre(2) ) .^ 2 ) );
sag = sqrt( R^2 - rho .^ 2) + s - R;

end
