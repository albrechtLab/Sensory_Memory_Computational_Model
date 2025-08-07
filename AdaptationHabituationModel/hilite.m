
function hilite(x0,y0,col,axlist,replace)
% function hilite(x,y,col,axlist,replace)
%       x is Nx2 matrix of [start end] x-pairs for rectangular highlight
%       y is Nx2 matrix of [start end] y-pairs for rectangular highlight
%           either x or y can be [] to span axis limits
%       color is 1x3 RGB color value (or [] = default)
%       axlist is list of axes ([] is current axis, or 'all' for all axes in current figure)
%       replace if true, clears prior highlighting (default = true)
%
% v1.1 2010-04-08
% Dirk Albrecht 

if ~exist('replace','var') || isempty(replace) replace = true; end
if ~exist('axlist','var') || isempty(axlist) axlist = gca; end
if ischar(axlist) axlist = findobj(gcf,'type','axes'); end
    
if ~exist('col','var') || isempty(col) col = [0.8 0.8 0.8]; end

if ~exist('x0','var') x0 = []; end 
if ~exist('y0','var') y0 = []; end 

% if ~isempty(isnan(x0)) && any(isnan(x0)) axlist = []; end  % no hilighting if NaN
% if ~isempty(isnan(y0)) && any(isnan(y0)) axlist = []; end  % no hilighting if NaN

if ~isempty(isnan(x0)) & any(isnan(x0)) axlist = []; end  % no hilighting if NaN
if ~isempty(isnan(y0)) & any(isnan(y0)) axlist = []; end  % no hilighting if NaN

for a = 1:length(axlist)
    ax = axlist(a);
    axlim = [get(ax,'XLim'), get(ax,'YLim')];
    
    if isempty(x0) x = axlim(1:2); else x = x0; end
    if isempty(y0) y = axlim(3:4); else y = y0; end

    % clear prior highlighing
    if replace
        h = findobj(ax,'Tag','hilite');
        if ~isempty(h) delete(h); end
    end
    
    axchildren = get(ax,'Children');

    h = [];
    for i = 1:size(x,1)
        for j = 1:size(y,1)
            h(i,j) = patch(x(i,[1 2 2 1 1]),y(j,[1 1 2 2 1]),col);
        end
    end

    set(h,'LineStyle','none','Parent',ax,'Tag','hilite');

    % order highlighting to back
    newhandles = reshape(h,[],1);
    set(ax,'Children',[axchildren; newhandles]);
end
    