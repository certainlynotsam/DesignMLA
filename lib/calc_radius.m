function R1 = calc_radius(f_mla, n_lens, type)
% Uses the lens maker's equation to determine the radius of curvature for a
% lens assuming a thin lens. 

switch type
    case 'plano'

    R1 = f_mla * ( n_lens - 1 ); % assumes R2 = 0

    otherwise
        warning('Unsupported lens type!')

end

end



