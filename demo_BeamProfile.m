%% A demo file for the BeamProfile class
%

%% Gaussian beams
% parse and check input
opt               = inputParser;
opt.StructExpand  = true;
opt.KeepUnmatched = false;
% number of pixels in x
opt.addParameter('nX', 1024, ...
    @(x) isnumeric(x) && isscalar(x));
% number of pixels in y
opt.addParameter('nY', 1360, ...
    @(x) isnumeric(x) && isscalar(x));
% pixel resolution in m/pix
opt.addParameter('pixres', 6.45e-6, ...
    @(x) isnumeric(x) && isscalar(x));
% beam waist at in m
opt.addParameter('w0', 200e-6, ...
    @(x) isnumeric(x) && isscalar(x));
% beam waist at in m
opt.addParameter('z', -1:0.1:1, ...
    @(x) isnumeric(x) && isscalar(x));
% energy in J
opt.addParameter('e', 0.2, ...
    @(x) isnumeric(x) && isscalar(x));
% lambda in m
opt.addParameter('lambda', 532e-9, ...
    @(x) isnumeric(x) && isscalar(x));
opt.parse();
opt = opt.Results;
% create fluence
x        = ((1:opt.nX)-(1+opt.nX)/2)*opt.pixres;
y        = ((1:opt.nY)-(1+opt.nY)/2)*opt.pixres;
[xg, yg] = ndgrid(x,y);
P   = @(z) exp(-2*((0.87*xg).^2+yg.^2)/(opt.w0*sqrt(1+(z*opt.lambda/(pi*opt.w0^2))^2))^2);
obj = BeamProfile(numel(opt.z));
for i = 1:numel(obj)
    obj(i).wavelength  = opt.lambda;
    obj(i).pixres      = opt.pixres;
    obj(i).cdata.cdata = P(opt.z(i));
    obj(i).energy      = opt.e;
    obj(i).z           = opt.z(i);
    
end
obj.iso11146;
obj.plotPropagation;
obj(1).plotProfile;