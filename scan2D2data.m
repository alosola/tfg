function [M_out] = scan2D2data(img_in,cmap,data_range,nan_marker)
% SCAN2D2DATA interpret pseudocolor plots, estimating the data.
%
% Calling:
% [M_out] = scan2D2data(img_in,cmap,data_range,nan_marker)
%
% Input
%   IMG_IN - rgb image containing pseudocolor plots
%   CMAP   - colour map of the pseudocolor plot, optional. If not
%            and given IMG_IN holds a colorbar it can be selected
%            interactively.
%   DATA_RANGE - [min max] to scale the data between.
%   NAN_MARKER - [r g b] color triplet used to mark nans.


%       Bjorn Gustavsson 1999-07-02
%	Copyright (c) 1999 by Bjorn Gustavsson

if nargin < 4 || isempty(nan_marker)
  nan_marker = [1,1,1];
end
nan_marker
if nargin <3
  [data_min] = input('Now give me the min data value:');
  [data_max] = input('Now give me the max data value:');
  data_range = [data_min data_max];
end
if nargin <= 2 || isempty(cmap)
  
  imagesc(img_in)
  
  set(gca,'position',[.05 .05 .9 .9])
  disp('Give me 2 diagonal corners containing the ''data region''')
  
  [x,y] = ginput(2);
  
  y = round([min(y) max(y)]);
  x = round([min(x) max(x)]);
  y(1) = max(1,y(1));
  x(1) = max(1,x(1));
  x(2) = min(x(2),size(img_in,2));
  y(2) = min(y(2),size(img_in,1));
  
  disp(size(img_in))
  disp(y)
  disp(x)
  
  M = img_in(y(1):y(2),x(1):x(2),:);
  
  if ( nargin == 1 ) || isempty(cmap)
    
    disp('Select the image region that contains the ''colormap''')
    
    [x,y] = ginput(2);
    
    y = round([min(y) max(y)]);
    x = round([min(x) max(x)]);
    
    cmap = (img_in(y(1):y(2),x(1):x(2),:));
    
    csize = size(cmap);
    
    if csize(1) > csize(2)
      cmap = median(double(cmap),2);
    else
    cmap = median(double(cmap),1);
    end
    cmap = squeeze(cmap);
    cmap = cmap/256;
    if max(cmap(:))>1
      cmap = cmap/max(cmap(:));
    end
    rgbplot(cmap),colormap(cmap),colorbar
    cm_proc = 'n';
    [cm_proc]=input('Do you need to flip the colormap around, y/n [n] :','s');
    if strcmp(cm_proc,'y')
      cmap = flipud(cmap);
    end
    cm_med = -1;
    [cm_med] = input('Do you want to smooth the colormap, 1 - (no smoothing) 8 (much) [1]: ');
    while ( cm_med > 0 )
      
      Cm_med = cm_med;
      cmap = medfilt1(cmap,Cm_med);
      rgbplot(cmap),colormap(cmap),colorbar
      [cm_med] = input('Smooth the colormap more/less, quit with -1: ');
      
    end
    
  end
    
    if any(exist('rgb2ind')==[2 3 5 6])
      M_out = rgb2ind(M,cmap,'nodither');
    else
      
      M1 = M(:,:,1);
      M2 = M(:,:,2);
      M_out = M(:,:,3);
      
      lmin = distmat1(double([M1(:),M2(:),M_out(:)]),...
                      [cmap(:,1)*256,cmap(:,2)*256,cmap(:,3)*256]);
      [M1,M2] = min(lmin');
      M_out(:) = M2(:);
    end
    %keyboard
    [data_min] = input('Now give me the min data value:');
    [data_max] = input('Now give me the max data value:');
    M_out = double(M_out);
    
    M_out = (data_max-data_min)/(length(cmap))*(M_out)+data_min;
    M_out = flipud(M_out);
    
elseif nargin >= 3
  
  if ndims(cmap) == 3
    
    % Assuming the colormap comes directly from image as NxMx3
    % 0-255. Here averaging over M, scaling to 0-1 and slightly
    % smoothing.
    cmap = filtfilt([1 1 1]/3,1,fliplr(squeeze(mean(cmap(:,:,:),2))))/256;
    
  end
  
  M1 = img_in(:,:,1);
  M2 = img_in(:,:,2);
  M_out = img_in(:,:,3);
  if min(size(nan_marker)==1)
    nan_indx = ( M1(:)==nan_marker(1) & ...
                 M2(:)==nan_marker(2) & ...
                 M_out(:)==nan_marker(3) );
  else % RGB-range
    nan_indx = ( nan_marker(1,1) < M1(:) & M1(:) < nan_marker(2,1) & ...
                 nan_marker(1,2) < M2(:) & M2(:) < nan_marker(2,2) & ...
                 nan_marker(1,3) < M_out(:) & M_out(:) < nan_marker(2,3) );
    disp(numel(nan_indx))
  end
  lmin = distmat1(double([M1(:),M2(:),M_out(:)]),...
                  [cmap(:,1)*256,cmap(:,2)*256,cmap(:,3)*256]);
  [M1,M2] = min(lmin');
  M_out(:) = M2(:);
  M_out = double(M_out);
  M_out(nan_indx) = nan;
  M_out = (max(data_range)-min(data_range))/(length(cmap))*(M_out)+min(data_range);
  M_out = flipud(M_out);
  
end
