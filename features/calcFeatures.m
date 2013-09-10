function [F,areas, Neigh] = calcFeatures(surface1, options)
  F = [];
  areas = [];
  if(~isfield(options,'method'))
    options.method = 'asi-hks';
  end
  if(~isfield(options,'hksscale'))
    options.hksscale = linspace(65,90,10);
    % options.hksscale = linspace(70,100,10);
    % options.hksscale = linspace(100,200,10);
  end
  if(~isfield(options,'laplacian'))

    options.hksscale = linspace(30,50,10);
    options.laplacian = 'graph';%'pcd'; %'cot'
  end

  if(~isfield(options,'h'))
    options.h = 1.5;
  end

  N1 = length(surface1.X(:));

  fprintf('%s\n',options.laplacian);
  switch(options.method)
    case 'x^p'

      %down sample
      if(size(Dist1,1) == 0)
        [sample1 Dist1] = farptsel(surface1.TRIV, surface1.X,surface1.Y,surface1.Z,N1);

      else
        sample1 = find(Dist1 == 0);
        b =  0:size(Dist1,1):size(Dist1,2)*size(Dist1,1)-1;
        sample1 = sample1 - b';
      end


      %create distance matrix  to contain only distances to sample
      F(:,[1:N1]) = Dist1(sample1,[1:N1]);

      opt = [];
      if(isfield(options,'center'))
        opt.center = 1;
      end
      if(isfield(options,'invert'))
        opt.invert= 1;
      end
      [W1, A1]  = mshlp_matrix(surface1,struct('dtype','cotangent')); %geodesic
      clear W1;clear W2;
      areas = [A1(sample1)];
    case 'hks-long'
      %calculate Laplacian

      [W1, A1]  = mshlp_matrix(surface1,struct('dtype','cotangent')); %geodesic


      %eigendecomposition
      Am1 = sparse([1:length(A1)], [1:length(A1)], A1);

      warning off;
      [evecs1 evals1] = eigs(W1, Am1, 300,-1e-5,struct('disp',0));
      evals1 = diag(evals1);
      warning on;
      evals1 = abs(evals1);
      L = 10;
      %   tn = [80 90 100];
      %tn = [ 16 32 48 64 72 ];
      %tn = linspace(60,90,L); %cutted
      %      tn = linspace(60,200,L);
      %            tn = [1000,1500];
      %           tn = [135.1 178.3 235.3 310.4];
      tn = [1351 1783 2353 3104]; % full
      %compute hks descripor (normalize by time also so that larger
      %times won't get less weight (too close to zero)
      % F = bsxfun(@times, (evecs1.^2)*exp(-bsxfun(@times, evals1(:), tn(:)')), (log(tn(:)').^2));
      %scale = (ones(size(evecs1)))*exp(-bsxfun(@times, evals1(:), tn(:)'));


      F =  (evecs1.^2)*exp(-bsxfun(@times, evals1(:), tn(:)'));
      %  F = F ./ scale;


      F  = normalize(F, 'l2', 2);

      %Give more weight to larger areas (vertex that belongs to
      %larger triangle will get more weight)

      %A1 = A1/(sum(A1));%+sum(A2));
      %A2 = A2/(sum(A2));%+sum(A2));
      areas = A1;
      %   A1 = A1/sum(A1);
      %         A2 = A1/sum(A2);
      %         F = bsxfun(@times, F, A1);
      %        M2 = bsxfun(@times, M2, 1+A2);     
    case 'hks'
      %calculate Laplacian



      %eigendecomposition


      warning off;
      switch (options.laplacian)
        case 'pcd'
          %[evecs1 evals1] = eigs(W1, Am1, 300,-1e-5,struct('disp',0));
          [W1,PcdNorm,Neigh,A1] = pcdlp_matrix([surface1.X(:),surface1.Y(:),surface1.Z(:)],2,struct('nn',10,'htype','psp','hs',1,'rho',1));
          %
          V = [surface1.X(:) surface1.Y(:) surface1.Z(:)]';
          N1 = length(surface1.X(:));
          for xi = 1:N1

            t = length(find(Neigh(xi,:) ~= 0))/2;%!t = size(ct,1)
            al = zeros(1,t);
            Amix = 0;
            for j = 1:t
              xj = Neigh(xi,j*2-1);
              xk = Neigh(xi,j*2);

              Xj = [V(:,xj)-V(:,xi) V(:,xk)-V(:,xi)];
              tarea = 0.5*sqrt(det(Xj'*Xj));
              al(j) = tarea;
              l0 = norm(V(:,xj) - V(:,xk));
              l1 = norm( V(:,xi) - V(:,xk));
              l2 = norm( V(:,xj) - V(:,xk));

              h0 = tarea / l0;
              h1 = tarea / l1;
              h2 = tarea / l2;

              cot0 = sqrt(abs(l1 * l1 - h2 * h2)) / h2;
              if( (l1 * l1 + l2 * l2 - l0 * l0) < 0)
                cot0= -cot0;
              end 
              cot1 = sqrt(abs(l2 * l2 - h0 * h0)) / h0;
              if( (l0 * l0 + l2 * l2 - l1 * l1) < 0)
                cot1 = -cot1;
              end
              cot2 = sqrt(abs(l0 * l0 - h1 * h1)) / h1;
              if( (l0 * l0 + l1 * l1 - l2 * l2) < 0 )
                cot2 = -cot2;
              end

              if( cot0 >= 0 && cot1 >= 0 && cot2 >= 0 )
                area0 = 1.0 / 8.0 * (l1 * l1 * cot1 + l2 * l2 * cot2);
                area1 = 1.0 / 8.0 * (l0 * l0 * cot0 + l2 * l2 * cot2);
                area2 = 1.0 / 8.0 * (l1 * l1 * cot1 + l0 * l0 * cot0);
              end
              if(cot0 < 0)
                area0 = tarea / 2;
                area1 = tarea / 4;
                area2 = tarea / 4;
              else  
                if(cot1 < 0)
                  area0 = tarea / 4;
                  area1 = tarea / 2;
                  area2 = tarea / 4;
                else

                  area0 = tarea / 4;
                  area1 = tarea / 4;
                  area2 = tarea / 2;

                end
              end 
              Amix = Amix + area0;
            end

            A1(xi) = sum(al)/3;
            %A1(xi) = Amix*2;
          end

          size(Neigh)
          Am1 = sparse([1:length(A1)], [1:length(A1)], A1);
          [evecs1 evals1] = eigs(W1,Am1, 300,-1e-5,struct('disp',0));
        case 'cot' 
          [W1, A1]  = mshlp_matrix(surface1,struct('dtype','cotangent')); %geodesic
          Neigh = [];
          %compute neigh
          nv = length(surface1.X(:));
          Neigh = sparse(nv,nv);
          for xi = 1:nv
            ct = getConnectedTris(surface1.TRIV, xi);
            ni = 1;
            j  = 1;
            for j = 1:length(ct)  
              tmp = find(surface1.TRIV(ct(j),:) ~= xi);
              xj = surface1.TRIV(ct(j),tmp(1));
              Neigh(xi,ni) = xj;
              ni = ni + 1;
              if(numel(tmp) > 1)
                xk = surface1.TRIV(ct(j),tmp(2));
                Neigh(xi,ni) = xk;
                ni = ni +1;
              end
            end
          end


          Am1 = sparse([1:length(A1)], [1:length(A1)], A1);
          [evecs1 evals1] = eigs(W1, Am1, 300,-1e-5,struct('disp',0)); %300
        case 'mesh'
          [W1, A1]  = mshlp_matrix(surface1, struct('dtype','euclidean', 'hs', 2, 'rho', 2, 'htype','ddr')); 

          Neigh = [];
          %compute neigh
         

          Am1 = sparse([1:length(A1)], [1:length(A1)], A1);
          [evecs1 evals1] = eigs(W1, Am1, 20,-1e-5,struct('disp',0)); %300
        case 'graph'
          tic
          [W1,PcdNorm,Neigh,A1] = pcdlp_matrix([surface1.X(:),surface1.Y(:),surface1.Z(:)],2,struct('nn',10,'htype','psp','hs',1,'rho',1));


          %perfect boundary
          %find edges
          %edge is edge when it belongs only to one triangle
          %             E = [surface1.TRIV(:,1) surface1.TRIV(:,2); 
          %             surface1.TRIV(:,1) surface1.TRIV(:,3);
          %             surface1.TRIV(:,2) surface1.TRIV(:,3)];
          %             E = sort(E,2);
          %             [E m n] = unique(E,'rows');
          %             c = accumarray(n,1);
          N = size(surface1.X,1);
          %             boundary = E( c == 1,:);
          %             boundary = unique(boundary); %% this is the vertex idxs of the edges
          %             %
          %
          %boundary = getBoundary(surface1);  
          boundary=[];
          %length(boundary)
          %change to 5
%          global decc;
          h = 2.7;
%          if exist('decc', 'var') && ~isempty(decc)
%            h = h * decc;
%          end
%          if N < 400
%            h = 100
%          end
          [W1 A1] = graphlp_matrix(surface1, 8, 6.25  ,boundary,A1);
%          [W1 A1_] = graphlp_matrixold(surface1, 8, 1.5 ,boundary);
          toc
          fprintf('Laplacian calculated')
          nv = length(surface1.X(:));
          Am1 = sparse([1:nv], [1:nv], A1);
          if N > 320
            neiges = 300;
          else
            neiges = N - 10; 
          end
          tic
          [evecs1 evals1] = eigs(W1, inv(Am1), neiges,-1e-5,struct('disp',0));
          evecs1 =  inv(Am1)*evecs1;% * evecs1;
          toc
          fprintf('eigen decomposition ended');
          %             [evecs1 evals1] = eigs(W1, Am1, 300,-1e-5,struct('disp',0));
        end
        evals1 = diag(evals1);
        warning on;
        evals1 = abs(real(evals1));
        evecs1 = real(evecs1);
        %L = 10;
        %   tn = [80 90 100];
        %tn = [ 16 32 48 64 72 ];
        %tn = linspace(65,90,L); %cutted
        tn = options.hksscale;
        %      tn = linspace(60,200,L);
        %            tn = [1000,1500];
        %           tn = [135.1 178.3 235.3 310.4];
        % tn = [1351 1783 2353 3104]; % full
        %compute hks descripor (normalize by time also so that larger
        %times won't get less weight (too close to zero)
        % F = bsxfun(@times, (evecs1.^2)*exp(-bsxfun(@times, evals1(:), tn(:)')), (log(tn(:)').^2));
        %scale = (ones(size(evecs1)))*exp(-bsxfun(@times, evals1(:), tn(:)'));


        F =  (evecs1.^2)*exp(-bsxfun(@times, evals1(:), tn(:)'));
        %  F = F ./ scale;


        F  = normalize(F, 'l2', 2);

        %Give more weight to larger areas (vertex that belongs to
        %larger triangle will get more weight)

        %A1 = A1/(sum(A1));%+sum(A2));
        %A2 = A2/(sum(A2));%+sum(A2));
        areas = A1;
        %   A1 = A1/sum(A1);
        %         A2 = A1/sum(A2);
        %         F = bsxfun(@times, F, A1);
        %        M2 = bsxfun(@times, M2, 1+A2);
      case 'fpfh'
        [W1,PcdNorm,Neigh,A1] = pcdlp_matrix([surface1.X(:),surface1.Y(:),surface1.Z(:)],2,struct('nn',10,'htype','psp','hs',1,'rho',1));

        [F Norm2] = pclfeature([surface1.X(:),surface1.Y(:),surface1.Z(:)]', 15);
        F = F';
        F  = normalize(F, 'l1', 2);

        areas = A1;%ones(size(A1));
      case 'pfh'
        [W1,PcdNorm,Neigh,A1] = pcdlp_matrix([surface1.X(:),surface1.Y(:),surface1.Z(:)],2,struct('nn',10,'htype','psp','hs',1,'rho',1));

        [tmp Norm2 F] = pclfeature([surface1.X(:),surface1.Y(:),surface1.Z(:)]', 10);
        F = F';
        F  = normalize(F, 'l1', 2);

        areas = A1;%ones(size(A1));
      case 'wks'
        [W1,PcdNorm,Neigh,A1] = pcdlp_matrix([surface1.X(:),surface1.Y(:),surface1.Z(:)],2,struct('nn',10,'htype','psp','hs',1,'rho',1));


        N = size(surface1.X,1);
        %boundary = getBoundary(surface1);
        boundary=[];
        %length(boundary)
        [W1 A1] = graphlp_matrix(surface1, 8, 1.5,boundary,A1);
        nv = length(surface1.X(:));
        Am1 = sparse([1:nv], [1:nv], A1);
        %W1 = W1 * inv(Am1);
        %W1 = (W1 + W1')/2;
        if N > 320
          neiges = 300
        else
          neiges = 200
        end
        tic
        [evecs1 evals1] = eigs(W1, inv(Am1), neiges,-1e-5,struct('disp',0));
        size(inv(Am1) * evecs1)
        evecs1 =  inv(Am1)*evecs1;% * evecs1;
        evals1 = diag(evals1);
        warning on;
        evals1 = abs(real(evals1));
        sum(imag(evecs1(:)))
        evecs1 = real(evecs1);

        evals1 = abs(evals1(:));  
        % calculate WKS
        emin = min(log(abs(evals1(:))))
        emax = max(log(abs(evals1(:))))
        s = 0.5*(emax - emin)/neiges
        emin = emin;%2*s;
        emax = emax;
        %           emin = -1000;
        %           emax = 0;
        tn = linspace(emin, emax, 1000);

        ce = sum( exp( -((bsxfun(@minus, log(evals1(:)) , tn(:)')).^2)./(2*s*s)));
        ce(ce == 0) = 1;

        ce = 1./ce;

        tmp = bsxfun(@minus, log(evals1(:)), tn(:)');
        size(tmp)
        size(evecs1)
        lognorm = (evecs1.^2)*exp( -(tmp.^2)./(2*s*s));
        F = lognorm * diag(ce);
        F  = normalize(F, 'l1', 2);
      case 'si-hks'
        %calculate Laplacian
        [W1, A1]  = mshlp_matrix(surface1,struct('dtype','cotangent')); %geodesic

        %eigendecomposition
        Am1 = sparse([1:length(A1)], [1:length(A1)], A1);

        warning off;
        [evecs1 evals1] = eigs(W1, Am1, 200,-1e-5,struct('disp',0));
        evals1 = diag(evals1);
        warning on;
        evals1 = abs(evals1);

        tau =linspace(-7,17,1024);
        tn = exp(tau);

        %compute hks descripor
        h1 = log( (evecs1.^2)*exp(-bsxfun(@times, evals1(:), tn(:)')));
        h1 = log( (evecs1.^2)*exp(-bsxfun(@times, evals1(:), exp(tau+1)))) - h1;
        F = abs(fft(h1')); 
        F =F(1:20,:);F = F';



        F  = normalize(F, 'l2', 2);


        %Give more weight to larger areas (vertex that belongs to
        %larger triangle will get more weight)

        %A1 = A1/(sum(A1));%+sum(A2));
        %A2 = A2/(sum(A2));%+sum(A2));
        areas = A1;
        %  A1 = A1/sum(A1(:)) ;
        %          A2 = A2/sum(A2(:)) ;
        %F = bsxfun(@times,F, A1);
        %           M2 = bsxfun(@times, M2, A2);


      case 'asi-hks'
        %calculate Laplacian
        [W1, A1]  = mshlp_matrix(surface1,struct('dtype','cotangent')); %geodesic

        %eigendecomposition
        Am1 = sparse([1:length(A1)], [1:length(A1)], A1);

        warning off;
        [evecs1 evals1] = eigs(W1, Am1, 300,-1e-5,struct('disp',0));
        evals1 = diag(evals1);
        warning on;
        evals1 = abs(evals1);

        %tau =linspace(4.1,4.6,1024); %- cutted
        tau =linspace(-7,17,1024); %-% full
        %            tau = linspace(1,10,512);
        tn = exp(tau);

        %compute hks descripor
        h1 = log( (evecs1.^2)*exp(-bsxfun(@times, evals1(:), tn(:)')));
        h1 =  ((evecs1.^2)*(bsxfun(@times, bsxfun(@times,exp(-bsxfun(@times, evals1(:), tn(:)')),-evals1(:)),tn(:)')))./h1;

        F = abs(fft(h1')); 
        F =F(1:20,:);F = F';


        %F = bsxfun(@times,F, A1);      
        F  = normalize(F, 'l2', 2);
        areas = A1;
      case 'asi-hks2'
        %calculate Laplacian
        [W1, A1]  = mshlp_matrix(surface1,struct('dtype','cotangent')); %geodesic

        %eigendecomposition
        Am1 = sparse([1:length(A1)], [1:length(A1)], A1);

        warning off;
        [evecs1 evals1] = eigs(W1, Am1, 200,-1e-5,struct('disp',0));
        evals1 = diag(evals1);
        warning on;
        evals1 = abs(evals1);

        eta =linspace(-20,20,1024);
        omega = exp(eta);

        %compute hks descripor
        h1 = (evecs1.^2)*bsxfun(@times,(-1./(( bsxfun(@plus, evals1,j*omega(:)').^2))),j*omega(:)');

        F = abs(ifft(real(h1'))); 
        F =F(1:20,:);F = F';



        F  = normalize(F, 'l2', 2);
        areas = A1;
      case 'vq' 
        %calculate Laplacian
        op = [];
        op.method='hks-long';
        [M1 areas] = calcFeatures(surface1, N1,Dist1,op);
        %%%%%%%%%

        vs = 48; % vocabulary size
        [idx1, vocab1] = kmeans(double(M1), vs, ...
          'maxiter', 300, ...
          'display', 1, ...
          'replicates', 1, ...
          'randstate', 0, ...
          'outlierfrac', 1e-3);

        % Compute sigma
        tree = ann('init', vocab1');
        [sig, mind] = ann('search', tree, M1', 1, 'eps', 1.1);
        sigma1 = median(mind);
        ann('deinit', tree);
        ann('close');

        F = {};
        F{1} = vocab1;
        F{2} = sigma1;

        %Give more weight to larger areas (vertex that belongs to
        %larger triangle will get more weight)

        %A1 = A1/(sum(A1));%+sum(A2));
        %A2 = A2/(sum(A2));%+sum(A2));
        % areas = [A1;A2];
        %M1 = bsxfun(@times, M1, A1');
        %M2 = bsxfun(@times, M2, A2');
        %F = F';


      otherwise
        error('unsupported method specified');
      end
    end



