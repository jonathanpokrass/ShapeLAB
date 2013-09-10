function [P,D,d] = farptsel(shape, N, varargin)
  TRIV = shape.TRIV;
  X = shape.X;
  Y = shape.Y;
  Z = shape.Z;

  Nv = length(X);

  defaultOpt.sample = round(rand*(Nv-1)) + 1; % Initialize with a random point.;
  defaultOpt.minRadius = 0;
  options = parseOpt(defaultOpt, varargin{:});

  f = fastmarchmex('init', int32(TRIV-1), double(X(:)), double(Y(:)), double(Z(:)));


  sample = options.sample;
  D      = repmat(Inf, [Nv N]);    % Distance maps.
  d      = repmat(Inf, [Nv 1]);
  P = zeros(N+1,1);
  if Nv ~= N
    P(1) = sample;
  else
    P(1) = 1;
  end
  idx = 1;
  for k = 1:N
    % Compute distance map from sample on the shape.
    u = repmat(Inf, [Nv 1]);
    u(P(k)) = 0;

    D(:,k) = fastmarchmex('march',f,double(u));

    d = min(d, D(:,k));
    
    if Nv ~= N
      [r, idx] = max(d);
    else
      idx = idx + 1;
    end

    P(k+1) = idx;
    if r < options.minRadius
        break;
    end
    %    fprintf('\b\b\b\b\b%3dpt',k)
  end
  P = P(1:min(k, N));
  D = D(:,1:min(k, N));
  %fprintf('\n',k)
  fastmarchmex('deinit', f);
end
