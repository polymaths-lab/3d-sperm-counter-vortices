function rotmat = vrrotvec2mat_custom(axisrot)


    axis = axisrot(1:3);
    theta = axisrot(4);
    axis = axis / norm(axis);
    x = axis(1); y = axis(2); z = axis(3);
    c = cos(theta);
    s = sin(theta);
    C = 1 - c;

    rotmat = [x*x*C + c,   x*y*C - z*s, x*z*C + y*s;
         y*x*C + z*s, y*y*C + c,   y*z*C - x*s;
         z*x*C - y*s, z*y*C + x*s, z*z*C + c];