function axisrot = vrrotvec_custom(dir_raw,dir_new)


    dir_raw = dir_raw / norm(dir_raw);
    dir_new = dir_new / norm(dir_new);
    axis = cross(dir_raw, dir_new);
    angle = acos(dot(dir_raw, dir_new));
    axis = axis / norm(axis);
    axisrot = [axis, angle];