% Function to draw a 3D rectangle
% TODO: make this into studied struct

function drawRect3D(position, dimensions, color)
    x = [0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0] * dimensions(1) + position(1);
    y = [0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1] * dimensions(2) + position(2);
    z = [0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1] * dimensions(3) + position(3);
    
    patch(x, y, z, color);
end