
Gaussian = @(x,sx)exp( -x.^2 / (2*sx^2));
% WeightMap = @(f0,sv)Gaussian( repmat(f0, [1 1 p])...
%     - repmat( reshape((0:p-1)/(p-1), [1 1 p]) , [n n 1]), sv );
prompt = 'sx?';
sx = input(prompt);
prompt = 'x?';
x = input(prompt);

g = Gaussian(x, sx);
