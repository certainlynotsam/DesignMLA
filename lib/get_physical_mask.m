function[thickness] = get_physical_mask(D_mla,D_cam,N_pixels_obj,N_pixels_mla,x_mla,y_mla,curvature)
% Compute sagitta for the 3-hex MLA. 

R = curvature; % radius of curvature
interval = D_cam; % powerphotonics xy-grid

% Generate MLA + mask
N_views = numel(x_mla); % total number of microlenses to be processed
[px, py] = meshgrid((-N_pixels_obj/2):(N_pixels_obj/2-1),(-N_pixels_obj/2):(N_pixels_obj/2-1)); % construct xy-grid of points equal to size of conj. BFP
u = zeros(size(px,1),size(px,2)); v = u; % initialize (intermediate) u and v matrices 
mla_pattern = zeros(N_pixels_obj,N_pixels_obj); % initialize matrix the size of the conj. BFP
total_inpolygon = zeros(size(px,1),size(px,2)); % initialize (intermediate) matrix for counting hexagons later

for i = 1 : N_views
    hexgon_shape = nsidedpoly(6,'Center',[x_mla(i) y_mla(i)],'SideLength',N_pixels_mla/sqrt(3)); % create a hexagon centred on a given u-lens coordinate
    points_on_hexagon_X(i,:) = (hexgon_shape.Vertices(:,1)-x_mla(i))*cos(pi/6)-(hexgon_shape.Vertices(:,2)-y_mla(i))*sin(pi/6)+x_mla(i); % rotate hexagon 30 degrees
    points_on_hexagon_Y(i,:) = (hexgon_shape.Vertices(:,2)-y_mla(i))*cos(pi/6)+(hexgon_shape.Vertices(:,1)-x_mla(i))*sin(pi/6)+y_mla(i);

    in_polygon = inpolygon(px, py, points_on_hexagon_X(i,:), points_on_hexagon_Y(i,:)); % boolean - choose points in px and py that lie within the hexagon
    u = u + in_polygon.*x_mla(i) / N_pixels_mla;
    v = v + in_polygon.*y_mla(i) / N_pixels_mla;
    lens_shape = calc_sag(R, N_pixels_obj, N_pixels_mla, interval, [x_mla(i) y_mla(i)]); % calculate sagittus at given coordinates
    lens_shape = lens_shape.*in_polygon; % mask the lens shape to only apply it within the hexagonal area
    mla_pattern = mla_pattern+lens_shape; % add the lens shape to the overall MLA pattern
    total_inpolygon = total_inpolygon+in_polygon; % add to count of used hexagons
end

thickness = mla_pattern;

%% Legacy 
% % Smoothing effect
% % c_matrix = ones(20,20)/400;
% % mask_pattern = conv2(mask_pattern,c_matrix);
% % new_range = 10:(size(mask_pattern,1)-10);
% % mask_pattern = mask_pattern(new_range,new_range);
% % thickness = thickness - mask_pattern - sag;
% thickness = mla_pattern;
% % thickness = (max(thickness,[],'all')-thickness).*total_inpolygon;
% 
% % thickness = thickness - not(total_inpolygon)*max(mask_pattern+sag,[],'all');

