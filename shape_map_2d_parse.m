function [shape_array,shape_polygon] = shape_map_2d_parse(dim,varargin)
        
%xdim and ydim values specify the sidelengths
%Will have to convert to scale-factor for SLM pattern updates

%REQUIRED INPUT: dimension (density of printing)
%OPTIONAL INPUTS: x/y area size, starting coordinates,shape to print
shape_polygon = 0;

defaultX = 5e-06;
defaultY = defaultX;
defaultCoords = [0,0];
defaultShape = 'rect';
defaultCirc = defaultX;

p = inputParser;
p.addParameter('xdim',defaultX);
p.addParameter('ydim',defaultY);
p.addParameter('dim',dim);
p.addParameter('coords',defaultCoords);
p.addParameter('shape',defaultShape);
p.addParameter('circle',defaultCirc);

defaultImage = 'fix_this_test2';
p.addOptional('image_name',defaultImage)
p.parse(varargin{:});

%TO DO: ADD INPUT SPECIFICATION REQUIREMENTS

shape_print = p.Results.shape;
co_input = [p.Results.coords(1),p.Results.coords(2)];

switch shape_print
    case 'rect'
        disp('Generating rectangle bitmap');
        [shape_array,shape_polygon] = generate_rectangle(co_input,p.Results.xdim,p.Results.ydim,p.Results.dim);
    case 'triangle'
        disp('Generating triangle bitmap');
        [shape_array,shape_polygon] = generate_triangle(co_input,p.Results.xdim,p.Results.ydim,p.Results.dim);
    case 'circle'
        disp('Generating circle bitmatp');
        [shape_array] = generate_circle(co_input,p.Results.xdim,p.Results.dim); %xdim=circumference
        shape_polygon = 1;
    case 'image'
        disp('Arbitary image conversion');
        parse_image = p.Results.image_name;
        [shape_array] = generate_bitmap_arb(co_input,p.Results.xdim,p.Results.ydim,p.Results.dim,parse_image);
        shape_polygon = 1;
  otherwise
        disp('Error.....');
end


end