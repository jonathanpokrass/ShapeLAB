% Convert argument pairs to a structure
function S = args2struct(varargin)

% No inputs, return empty structure
if isempty(varargin), S = struct(); return; end

if length(varargin) == 1 && isempty(varargin{1}), S = struct(); return; end

% Need pairs of inputs
if mod(length(varargin),2)==1
   error('number of arguments must be even');
end

% Odd elements of varargin are fields, even ones are values
% Store all field names in lower case
for k = 1:2:length(varargin)
   S.(varargin{k}) = varargin{k+1};
end

