% 
% Finds the point to point correspondence between shape1 and shape2
% from C, so that
% C * basis1 ~ basis2(:, shape1toshape2) 
% C' * basis2 ~ basis1(:, shape2toshape1)
% Author Jonathan Pokrass
function [shape1toshape2, shape2toshape1, refinedC] = ...
         calcP2PFromC(shape1, shape2, C, basis1, basis2, varargin)
defaultOpt.useGroundTruth = 1;
defaultOpt.useSymmetricToo = 1;
defaultOpt.numRefinements = 10;
defaultOpt.debug = 0;
defaultOpt.drawSphere = 0;

opt = parseOpt(defaultOpt, varargin{:});

C = C';
refinedC = C;
prevC = C;
timeit.ICP = 0;
b2 = basis2';
b1 = basis1';
prevVal = Inf;
searchIndexParams = struct();  
shape1toshape2 = [];
for icpIter = 0:opt.numRefinements
  if opt.debug
    figure(1)
    subplot(1,2,1);
    imagesc(abs(C))
    subplot(1,2,2);
    imagesc(abs(refinedC))
    title('refined C');
  end
  
  b2Perm = refinedC * b1;
  b1Perm = refinedC'* b2;
  if icpIter == opt.numRefinements
    %searchIndexParams = struct('algorithm', 'linear');
  else
    searchIndexParams = struct();  
  end
  
  tic;
  %
  prevShape1toshape2 = shape1toshape2;
  shape1toshape2 = flann_search(b2, b2Perm, 1, searchIndexParams);
  %shape2toshape1 = flann_search(refinedC * basis1', b2, 1, searchIndexParams);
  B = b2(:, shape1toshape2);
  A = b1;
  newVal = norm(refinedC * A - B,'fro');
  if prevVal < newVal
    printlog('   Refine worsened the results restoring prevC');
    refinedC = prevC;
    shape1toshape2 = prevShape1toshape2;
    b2Perm = refinedC * b1;
    b1Perm = refinedC'* b2;
    break;
  end
  prevVal = newVal;
  
  flannTime = toc;
  timeit.ICP = timeit.ICP + flannTime;
  printlog('  flann time = %f secs', flannTime);
  if icpIter < opt.numRefinements
    printlog('  Preforming refinement (%d/%d)', icpIter + 1, opt.numRefinements);
    
    tic;
    %find refinedC = argmin ||refinedC * basis1 - basis2(:, shape1toshape2result)|| 
    % s.t. refinedC * refinedC' = I
    
    [U,S,V] = svd(B*A');
    U(:,end)=U(:,end)*det(U*V');
    newCR=U*V';
    timeit.ICP = timeit.ICP + toc;
    
    
    crdiff = sum(sum((abs(newCR - refinedC))));
    printlog('    crdiff = %f', crdiff);
    
    if crdiff < 0.01 
      printlog('  Refinments converged stopping refinments')
      break;
    end
    prevC = refinedC;
    refinedC = newCR;
  end
end
shape2toshape1 = flann_search(basis1', b1Perm, 1, searchIndexParams);

refinedC = refinedC';
%shape2toshape1result = flann_search(CR * basis1', b2, 1, searchIndexParams);
%flann_free_index(b2_index);

if opt.debug
  %Draw point-wise correspondance
  figure(10)
  clf
  subplot(1, 4, 1);
  
  shape1_ = rotate_shape(shape1, rotation_matrix(0,0,0*pi/180), [0 0 0]);
  shape2_ = rotate_shape(shape2, rotation_matrix(0,0,0*pi/180), [0 0 0]);

  d = shape1.X(:) + max(abs(shape1.X(:))) + 1;
  d2  = shape2.X(:) + max(abs(shape2.X(:))) + 1;

  nv1 = length(shape1toshape2);
  colormap(jet(1000));
  drawShape(shape1, 'ptColor', d([1:nv1]));
  title('     ');

  subplot(1, 4, 3);
  drawShape(shape2_, 'ptColor', d(shape2toshape1));
  title('result');
  colormap(jet(50));

  for j = 20:100:20
    points2idxs = farptsel(shape2, j);
    points1 = []; points2 = [];
    points2.X = shape2_.X(points2idxs);
    points2.Y = shape2_.Y(points2idxs);
    points2.Z = shape2_.Z(points2idxs);


    points1.X = shape1_.X(shape2toshape1(points2idxs));
    points1.Y = shape1_.Y(shape2toshape1(points2idxs));
    points1.Z = shape1_.Z(shape2toshape1(points2idxs));

    figure;
    colormap(jet(200));
    plotMatchingPoints(shape2_, shape1_, points2idxs, shape2toshape1(points2idxs),...
                   'ptColor1', d2, 'ptColor2', d2(shape1toshape2), 'drawSphere', opt.drawSphere);
    figure;
    plotMatchingPoints(shape1_, shape2_, shape2toshape1(points2idxs),points2idxs,...
                   'ptColor1', d, 'ptColor2', d(shape2toshape1), 'drawSphere', opt.drawSphere);
end
drawnow;
end
end