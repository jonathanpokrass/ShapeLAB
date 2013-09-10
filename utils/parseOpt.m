% Construct options structure with a template
% S can be varargin or struct
% D are the default values struct
function T = parseOpt(D,varargin)
if length(varargin) == 1 && isstruct(varargin{1}) 
	S = varargin{1};
else
	S = args2struct(varargin{:});
end
T     = D;             % copy the template
if isempty(S)
	return;
end
% Check arguments, must have two structures
if ~(isstruct(S) && isstruct(D))
   error('input arguments must be structures');
end
   
fname = fields(S);     % make a list of field names

% Loop over all fields in the template, copy matching values from S
for k = 1:length(fname)
   % Process matching field names in S
   if isfield(D,fname{k})
      % Is this a substructure ?
      if isstruct(T.(fname{k})) && isstruct(S.(fname{k}))
         % Recursively process the substructure
         T.(fname{k}) = parseOpt(T.(fname{k}),S.(fname{k}));
      % Not a substructure, copy field value from S
      else T.(fname{k}) = S.(fname{k});
      end
   end
end
end


