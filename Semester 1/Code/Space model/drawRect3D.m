% Function to draw a 3D rectangle
% TODO: make this into studied struct


function draw3DRectangle(position, dimensions, color)
    % Define the vertices of the rectangle
    x = dimensions(1);
    y = dimensions(2);
    z = dimensions(3);

    x0 = position(1);

    x0plus = x0 + x/2;
    x0minus = x0 - x/2;

    vert = ;%TODO - I need to better understand how patch works
    fac = ;%TODO

    patch =  ('Vertices',vert,'Faces',fac)   
end