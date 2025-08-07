function h = stateplot3(mat,sem,xval,sep,stair,codelist,dotted,cmap)

[numc,numx] = size(mat);

if nargin < 8
    cmap = [1 1 1; 0.9 0.9 0.9; .7 .7 .7; .3 .3 .3; 0 0 0; .6 0 0; 1 .2 .2; 0 0 1; 0 0 0.5];
end
if nargin < 7 dotted = 1; end
if nargin < 6 codelist = 1:numc; end
if nargin < 5 stair = 1; end
if nargin < 4 sep = 0; end
if nargin < 3 xval = 1:numx; end
if nargin < 2 sem = NaN; end

if isnan(sem) sem = zeros(size(mat)); end
    
sepsnap = 0.1;
gridsep = 0.1;

delx = mean(diff(xval));

mat2 = [zeros(1,numx); flipud(mat)];
sem2 = [zeros(1,numx); flipud(sem)];
codelist2 = sort(numc - codelist + 1);

base = ones(numc,1)*(max(abs(sep),sepsnap));
reset(gca);

for i = codelist2
    xpos = xval;
    
    if sep
        ydn = repmat(base(i),1,numx);
        yup = ydn + mat2(i+1,:);
        if length(sep) > 1
            base(i+1) = base(i) + sep(i);
        else            
            if sep < 0
                base(i+1) = max(max(yup))-sep;
            else
                base(i+1) = ceil((max(max(yup))+sep)/sepsnap)*sepsnap;
            end    
        end
        maxy = base(i+1);
    else
        unused = setdiff(1:numc,codelist2+1);
        mat2(unused,:)=0;
        ydn = sum(mat2(1:i,:),1);
        yup = sum(mat2(1:i+1,:),1);
    end
    
    if stair
        yup = reshape(repmat(yup,2,1),1,[]);
        ydn = reshape(repmat(ydn,2,1),1,[]);
        xpos = reshape(repmat(xval,2,1)+repmat([-1;1]*delx/2,1,numx),1,[]);
    end
    if sep
        if stair
            jbfill(xpos,yup,ydn,cmap(numc-i+2,:),cmap(numc-i+2,:),i>1,1);
        else
            jbfill(xpos,yup-sem2(i+1,:),yup+sem2(i+1,:),cmap(numc-i+2,:)/2 + 0.5,cmap(numc-i+2,:)/2 + 0.5,i>1,1);
        end
        hold on; 
        plot(xpos,yup,'Color',cmap(numc-i+2,:)/2);
    else
        jbfill(xpos,yup,ydn,cmap(numc-i+2,:),-1,i>1,1);
    end

end
axis tight

set(gca,'XLim',[min(xval) max(xval)] + (stair>0)*[-1 1]*delx/2);
if sep
    labels = ['Unk';'Fwd';'Pau';'Rev';'OmR';'OmF';'60F';'60R'];
    idx = sort(numc-codelist2);
    set(gca,'YTick',base(codelist2),'YTickLabel',flipud(labels(idx+1,:)));
    set(gca,'YLim',[0 maxy]);
else
    set(gca,'YLim',[0 1]);
end
hold on;

if sep
    ylist = [];
    for i = codelist2
        ylist = [ylist, base(i)+gridsep:gridsep:(base(i+1)-0.001)];
    end
    if dotted
        for i = 1:length(ylist)
            line(get(gca,'XLim'),[1 1]*ylist(i),'Color',[1 1 1],'LineStyle',':','LineWidth',0.2);
        end
    end
end
    

%% 
function[fillhandle,msg]=jbfill(xpoints,upper,lower,color,edge,add,transparency,stair)

if nargin<8;stair=0;end %default is to have smooth curve, not stairs
if nargin<7;transparency=.5;end %default is to have a transparency of .5
if nargin<6;add=1;end     %default is to add to current plot
if nargin<5;edge='k';end  %dfault edge color is black
if nargin<4;color='b';end %default color is blue

upper(find(isnan(upper)))=0;
lower(find(isnan(lower)))=0;
if edge == -1;
    edge = 'k';
    noedge = 1;
else
    noedge = 0;
end
if stair
    if diff(size(xpoints))<1 xpoints = xpoints'; end
    if diff(size(upper))<1 upper = upper'; end
    if diff(size(lower))<1 lower = lower'; end
    dx = mean(diff(xpoints));
    xpoints = reshape(addortho([-dx/2; dx/2],xpoints),[],1)';
    upper = reshape(repmat(upper,2,1),[],1)';
    lower = reshape(repmat(lower,2,1),[],1)';
end

if length(upper)==length(lower) && length(lower)==length(xpoints)
    msg='';
    filled=[upper,fliplr(lower)];
    xpoints=[xpoints,fliplr(xpoints)];
    if add
        hold on
    end
    fillhandle=fill(xpoints,filled,color);
    set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency,'EdgeAlpha',transparency);%set edge color
    
    if noedge
        set(fillhandle,'LineStyle','none');
    end
    
    if add
        hold off
    end
    set(gca,'Layer','top');
else
    msg='Error: Must use the same number of points in each vector';
end
