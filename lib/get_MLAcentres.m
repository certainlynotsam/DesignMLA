function[x,y]=get_MLAcentres(N_mla,N_pixels_mla)
% Calculate the coordinate for centre of each microlens
% assume seamless tiled microlens
% x-x positions of microlens
% y-y positions of microlens
% N_mla-number of microlens in one dimension
% size_plane-number of pixels in MLA plane, about 468 pixel for 3x3

if N_mla >= 2 && N_mla < 3
    R = N_pixels_mla / sqrt(3);
    shortR = R / 2 * sqrt(3);
    x = [0 -shortR shortR -2*shortR 2*shortR -shortR shortR -3*shortR 3*shortR 0 -2*shortR 2*shortR];
    y = [R -R/2 -R/2 R R 2.5*R 2.5*R -R/2 -R/2 -2*R -2*R -2*R];

elseif N_mla >= 3 && N_mla < 4
    R = N_pixels_mla / sqrt(3);
    shortR = R / 2*sqrt(3);
    x_mla = [0 0 0 1.5*R -1.5*R 1.5*R -1.5*R];
    y_mla = [0 -2*shortR 2*shortR shortR shortR -shortR -shortR];
    x_outside = [1.5*R -1.5*R 1.5*R -1.5*R 3*R 3*R 3*R -3*R -3*R -3*R];
    y_outside = [3*shortR 3*shortR -3*shortR -3*shortR 0 2*shortR -2*shortR 0 2*shortR -2*shortR];
    x = [x_mla,x_outside];
    y = [y_mla,y_outside];
end

end


        
   
    
