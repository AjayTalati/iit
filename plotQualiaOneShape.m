function [ax, height, extra_plots] = plotQualiaOneShape(whole_p, whole_f, cut_p, cut_f,nWholeConcepts, whole_purviews, part_purviews,...
                                            highlight_indices, parent_panel, dim_option, all_phi, unconstrained)
                                        
% BASED ON GPLOTMATRIX
Fancy = 1;
if Fancy == 1
    textCol = [1,1,1];
else
    textCol = [0,0,0];
end

view_option = '3D';
x = [whole_p' whole_f'; cut_p' cut_f'];


num_dims = min(size(x,2),8);
num_nodes = log2(size(x,2)/2);

% princomp(x)

if strcmp(dim_option,'Variance')
    concept_value = var(x);
elseif strcmp(dim_option,'Mode')
%     concept_value = max(x);
    concept_value = sum(x,1);
end

[ignore_var state_ordering_temp] = sort(concept_value,'descend');
state_ordering_f = state_ordering_temp(state_ordering_temp > 2^num_nodes/2);
state_ordering_p = state_ordering_temp(state_ordering_temp <= 2^num_nodes/2);
%choose the two future states with max variance and the past state with max
%variance
state_ordering = [state_ordering_f(1:2) state_ordering_p(1)];

% This works also for the future state > 8!
state_label = {};
for s = 1:length(state_ordering)
    state_label = [state_label; strcat(int2str(index2state(state_ordering(s), 2*ones(1,num_nodes)))')];
end

for s = setdiff(state_ordering_temp,state_ordering)  
 state_label = [state_label; strcat(int2str(index2state(s, 2*ones(1,num_nodes)))')];
end

nPartsConcepts = size(x,1) - nWholeConcepts;

whole = x(1:nWholeConcepts,:);
part = x(nWholeConcepts+1:end,:);

w_highlight_indices = highlight_indices(highlight_indices <= nWholeConcepts);
p_highlight_indices = highlight_indices(highlight_indices > nWholeConcepts) - nWholeConcepts;

w_nonhighlight_indices = setdiff(1:nWholeConcepts,w_highlight_indices);
p_nonhighlight_indices = setdiff(1:nPartsConcepts,p_highlight_indices);


if Fancy == 1
    whole_labels =  cell(nWholeConcepts,1);
    for i = 1:nWholeConcepts
        whole_labels{i} = char(whole_purviews{i}-1+'A');
    end
end

whole_selected_labels = cell(length(w_highlight_indices),1);
for i = 1:length(w_highlight_indices)
    whole_selected_labels{i} = mod_mat2str(whole_purviews{w_highlight_indices(i)});
end

whole_nonselected_labels = cell(length(w_nonhighlight_indices),1);
for i = 1:length(w_nonhighlight_indices)
    whole_nonselected_labels{i} = mod_mat2str(whole_purviews{w_nonhighlight_indices(i)});
end

part_selected_labels = cell(length(p_highlight_indices),1);
for i = 1:length(p_highlight_indices)
    part_selected_labels{i} = mod_mat2str(part_purviews{p_highlight_indices(i)});
end

part_nonselected_labels = cell(length(p_highlight_indices),1);
for i = 1:length(p_nonhighlight_indices)
    part_nonselected_labels{i} = mod_mat2str(part_purviews{p_nonhighlight_indices(i)});
end

dims = floor(size(x,2)/2);

rows = 7; 

pos = get(parent_panel,'Position');
width = pos(3)/rows;
height = width;
space = .04; % 2 percent space between axes

ax = cell(nchoosek(num_dims,2)+1,1); % all pairs of dims plus the 3D plot
ax_index = 1;


row_pos = ceil(rows/8) - .5;
col_pos = row_pos;
size_scale = 24*row_pos;
horiz_offset = .25;


% axPos = [row_pos*width+space-horiz_offset col_pos*height+space ...
%         width*(1-space)*size_scale height*(1-space)*size_scale];
axPos = [0,0.1,1,1]
axes3D = axes('Position',axPos, 'visible', 'on', 'Box','on','Parent',parent_panel,'DrawMode','fast');

ax{ax_index} = axes3D;

plot3(ax{ax_index},0,0,0,'+k');
hold on

%  x = r*sin(az).*cos(el);
%  y = -r*cos(az).*cos(el);
%  z = r*sin(el);
% phi1 = acos(-1 + 2*rand(n,1));
% th1 = 2*pi*rand(n,1);
% % Convert to cart.
% x = r1.*sin(phi1).*sin(th1) + X;
% y = r1.*sin(phi1).*cos(th1) + Y;
% z = r1.*cos(phi1) + Z;
% past other axes
if Fancy == 1
    for dd = 1:dims-1
        r1 = 1;
        color = [0.2,0.2,0.8];
        phi1 = acos(-1 + 2*rand);
        th1 = pi/2 + pi/2*rand;
        % Convert to cart.
        xS = r1.*sin(phi1).*sin(th1);
        yS = r1.*sin(phi1).*cos(th1);
        zS = r1.*cos(phi1);
        plot3(ax{ax_index},[0 xS],[0 yS],[0 zS],'--', 'Color', color)
        plot3(ax{ax_index},xS,yS,zS,'^', 'Color', color)
        text(xS-0.03,yS-0.03,zS,state_label(3+dd), 'Color', color) %the first three states are selected
        hold on
    end

    % future other axes
    for dd = 1:dims-2
        color = [0,0.5,0];
        phi1 = acos(-1 + 2*rand);
        th1 = 3/2*pi + pi/2*rand;
        % Convert to cart.
        xS = r1.*sin(phi1).*sin(th1);
        yS = r1.*sin(phi1).*cos(th1);
        zS = r1.*cos(phi1);
        plot3(ax{ax_index},[0 xS],[0 yS],[0 zS],'--', 'Color', color)
        plot3(ax{ax_index},xS,yS,zS,'^', 'Color', color)
        text(xS-0.03,yS-0.03,zS,state_label(3+dd), 'Color', color) %the first three states are selected
        hold on
    end
end

plot3(ax{ax_index}, [0,1], [0,0],[0,0], '-g', 'LineWidth', 2)
plot3(ax{ax_index}, [0,0], [0,1],[0,0], '-g', 'LineWidth', 2)
plot3(ax{ax_index}, [0,0], [0,0],[0,1], '-b', 'LineWidth', 2)

plot3(ax{ax_index}, [1], [0],[0], '>g')
text(1+0.07,0+0.07,0,state_label(1), 'Color', textCol)

plot3(ax{ax_index}, [0], [1],[0], '<g')
text(0+0.07,1+0.07,0,state_label(2), 'Color', textCol)

plot3(ax{ax_index}, [0], [0],[1], '^b')
text(0+0.07,0+0.07,1,state_label(3), 'Color', textCol)

% % past other axes
% for dd = 1:dims-1
%     ax_CO = pi + dd*(pi/2)/(dims-1);
%     el = dd*(pi/2)/(dims-1);
%     color = [0,0,1];
%     plot3(ax{ax_index},[0,-cos(ax_CO)*cos(el)],[0,sin(ax_CO)*cos(el)],[0,cos(el)],'--', 'Color', color)
%     plot3(ax{ax_index},-cos(ax_CO),sin(ax_CO),0,'^', 'Color', color)
%     text(-cos(ax_CO)-0.03,sin(ax_CO)-0.03,0,state_label(3+dd), 'Color', textCol) %the first three states are selected
%     hold on
%     
% end
% 
% % future other axes
% for dd = 1:dims-2
%     ax_CO = pi/2 + dd*(pi/2)/(dims-1);
%     color = [0,1,0];
%     plot3(ax{ax_index},[0, cos(ax_CO)],[0,sin(ax_CO)],[0,0],'--', 'Color', color)
%     plot3(ax{ax_index},cos(ax_CO),sin(ax_CO),0,'^', 'Color', color)
%     text(cos(ax_CO)-0.03,sin(ax_CO)-0.03,0,state_label(3+dims-1+dd), 'Color', textCol)
%     hold on
% end

%Unconstrained 
scatter3(ax{ax_index},unconstrained(state_ordering(1)),unconstrained(state_ordering(2)),...
    unconstrained(state_ordering(3)),'Marker','x','MarkerEdgeColor','k','SizeData',100,'Clipping','on')
hold on

for l = 1:nWholeConcepts
    if all_phi(2,l)==0
        color = [.5 .5 .5];
    else
        color = [1, 0, 0];
    end   
    scatter3(ax{ax_index},part(l,state_ordering(1)),part(l,state_ordering(2)),...
        part(l,state_ordering(3)),'Marker','d','MarkerEdgeColor',color,'SizeData',75*all_phi(2,l)*4,'Clipping','on')
    hold on
end
for l = 1:nWholeConcepts
    scatter3(ax{ax_index},whole(l,state_ordering(1)),whole(l,state_ordering(2)),...
        whole(l,state_ordering(3)),'Marker','p','MarkerFaceColor','y','MarkerEdgeColor',[0.5,0.5,0],'SizeData',250*all_phi(1,l)*4,'Clipping','on')
end
hold on

scatter3(ax{ax_index},whole(w_highlight_indices,state_ordering(1)),whole(w_highlight_indices,state_ordering(2)),...
    whole(w_highlight_indices,state_ordering(3)),'Marker','o','MarkerEdgeColor','g','SizeData',100,'Clipping','on')
hold on

scatter3(ax{ax_index},part(p_highlight_indices,state_ordering(1)),part(p_highlight_indices,state_ordering(2)),...
    part(p_highlight_indices,state_ordering(3)),'Marker','o','MarkerEdgeColor','g','SizeData',100,'Clipping','on')
hold on

xlabel(ax{ax_index},dec2bin(state_ordering(1)-1,num_nodes))
ylabel(ax{ax_index},dec2bin(state_ordering(2)-1,num_nodes))
zlabel(ax{ax_index},dec2bin(state_ordering(3)-1,num_nodes))


N = num_nodes;
maxEnt = [unconstrained(state_ordering)];
data = cell(2^N-1,3);
for s = 1:length(whole_purviews)
    ind = concept2index(whole_purviews{s},[1:N]);  
    data(ind,:) = num2cell(x(s,state_ordering));
end    
% make connections according to higher order concepts
for i = 1:N
        line(cell2mat([maxEnt(1,1), data(i,1)]), cell2mat([maxEnt(1,2), data(i,2)]), cell2mat([maxEnt(1,3), data(i,3)]), 'Color', [.7 .7 1] )       
        for j = setdiff(1:N, i)
            k = concept2index(sort([i,j]),[1:N]);
            %if numel(data([i k], :))
                line(cell2mat([data(i,1) data(k,1)]), cell2mat([data(i,2) data(k,2)]), cell2mat([data(i,3) data(k,3)]), 'Color', [.7 .7 1])
                %line([data(i,1) data(k,1)], [data(i,2) data(k,2)], [data(i,3) data(k,3)], 'Color', [.7 .7 1] )
            for m = setdiff(1:N, [i,j])
                n = concept2index(sort([i,j,m]),[1:N]);
                line(cell2mat([data(k,1) data(n,1)]), cell2mat([data(k,2) data(n,2)]), cell2mat([data(k,3) data(n,3)]), 'Color', [.7 .7 1] )
                %line([data(k,1) data(n,1)], [data(k,2) data(n,2)], [data(k,3) data(n,3)], 'Color', [.7 .7 1] )
                for o = setdiff(1:N, [i,j,m])
                    p = concept2index(sort([i,j,m,o]),[1:N]);
                    line(cell2mat([data(n,1) data(p,1)]), cell2mat([data(n,2) data(p,2)]), cell2mat([data(n,3) data(p,3)]), 'Color', [.7 .7 1] )
                    %line([data(n,1) data(p,1)], [data(n,2) data(p,2)], [data(n,3) data(p,3)], 'Color', [.7 .7 1] )
                end
            end
        end    
end


if Fancy == 1
    p1 = patch([-1.75 0 1.75 0],[0 -1.75 0 1.75], [-1.2 -1.2 -1.2 -1.2],'k');
    p2 = patch([-1.75 0 0 -1.75],[0 -1.75 -1.75 0], [1.2 1.2 -1.2 -1.2],'k');
    p3 = patch([-1.75 0 0 -1.75],[0 1.75 1.75 0], [1.2 1.2 -1.2 -1.2],'k');


    topcolor = [0 0 1];
    bottomcolor = [0 0 0];
    cdata(1,1,:) = bottomcolor;
    cdata(1,2,:) = topcolor;
    cdata(1,3,:) = [.5 .5 1];
    cdata(1,4,:) = topcolor;

    set(p1,'CData',cdata, ...
        'FaceColor','interp', ...
        'EdgeColor','none', ...
        'HitTest','off');


    topcolor = [0 0 1];
    bottomcolor = [0 0 0];
    cdata(1,1,:) = bottomcolor;
    cdata(1,2,:) = topcolor;
    cdata(1,3,:) = topcolor;
    cdata(1,4,:) = bottomcolor;

    set(p2,'CData',cdata, ...
        'FaceColor','interp', ...
        'EdgeColor','none', ...
        'HitTest','off');


    topcolor = [0.1 0.1 1];
    bottomcolor = [0 0 0.3];
    cdata(1,1,:) = bottomcolor;
    cdata(1,2,:) = topcolor;
    cdata(1,3,:) = topcolor;
    cdata(1,4,:) = [0 0 0.1];

    set(p3,'CData',cdata, ...
        'FaceColor','interp', ...
        'EdgeColor','none', ...
        'HitTest','off');
    
    text(whole(:,state_ordering(1)),whole(:,state_ordering(2)),whole(:,state_ordering(3)),whole_labels, 'Color', textCol)
end

set(ax{ax_index},'xlimmode','manual','ylimmode','manual',...
        'xlim',[-1.1 1.1],'ylim',[-1.1 1.2],'zlim',[-1.1 1.1],...
        'CameraPosition',[15,-1,5] ,'CameraViewAngleMode','manual')
% set(ax{ax_index},'xlimmode','manual','ylimmode','manual',...
%             'xlim',[-1 1],'ylim',[-1 1.0],'zlim',[-.25 1.25],...
%             'CameraViewAngleMode','manual','Clipping','on')

%     text(part(p_nonhighlight_indices,state_ordering(1)),part(p_nonhighlight_indices,state_ordering(2)),...
%         part(p_nonhighlight_indices,state_ordering(3)),part_nonselected_labels)
text(part(p_highlight_indices,state_ordering(1))+.03,part(p_highlight_indices,state_ordering(2)),...
    part(p_highlight_indices,state_ordering(3)),part_selected_labels, 'Color', textCol)
%     text(whole(w_nonhighlight_indices,state_ordering(1)),whole(w_nonhighlight_indices,state_ordering(2)),...
%         whole(w_nonhighlight_indices,state_ordering(3)),whole_nonselected_labels)
text(whole(w_highlight_indices,state_ordering(1)),whole(w_highlight_indices,state_ordering(2)),...
    whole(w_highlight_indices,state_ordering(3)),whole_selected_labels, 'Color', textCol)

end


% -----------------------------
function datatipTxt = gplotmatrixDatatipCallback(obj,evt)

target = get(evt,'Target');
ind = get(evt,'DataIndex');
pos = get(evt,'Position');

dtcallbackdata = getappdata(target,'dtcallbackdata');
[BigAx,gnum,row,col] = dtcallbackdata{:};

ginds = getappdata(BigAx,'ginds');
xnam = getappdata(BigAx,'xnam');
ynam = getappdata(BigAx,'ynam');
xdat = getappdata(BigAx,'x');
ydat = getappdata(BigAx,'y');
XvsX = getappdata(BigAx,'XvsX');
gn = getappdata(BigAx,'gn');

gind = ginds{gnum};
obsind = gind(ind);

xvals = xdat(obsind,:);
yvals = ydat(obsind,:);

x = xvals(col);
y = yvals(row);

if x~=pos(1) || y~=pos(2)
    % Something is inconsistent, display default datatip.
    datatipTxt = {sprintf('X: %s',num2str(pos(1))),sprintf('Y: %s',num2str(pos(2)))};
else
    if isempty(xnam)
        xnam = cell(size(xdat,2),1);
        for i = 1:size(xdat,2)
            xnam{i} = sprintf('xvar%s',num2str(i));
        end
    end
    if isempty(ynam)
        ynam = cell(size(ydat,2),1);
        for i = 1:size(ydat,2)
            ynam{i} = sprintf('yvar%s',num2str(i));
        end
    end

    % Generate datatip text.
    datatipTxt = {
        [xnam{col},': ',num2str(x)],...
        [ynam{row},': ',num2str(y)],...
        '',...
        sprintf('Observation: %s',num2str(obsind)),...
        };

    if ~isempty(gn)
        datatipTxt{end+1} = ['Group: ',gn{gnum}];
    end
    datatipTxt{end+1} = '';

    xnamTxt = cell(length(xvals),1);
    for i=1:length(xvals)
        xnamTxt{i} = [xnam{i} ': ' num2str(xvals(i))];
    end
    datatipTxt = {datatipTxt{:}, xnamTxt{:}};
    
    if ~XvsX
        ynamTxt = cell(length(yvals),1);
        for i=1:length(yvals)
            ynamTxt{i} = [ynam{i} ': ' num2str(yvals(i))];
        end
        datatipTxt = {datatipTxt{:}, ynamTxt{:}};
    end

end
end

function [ogroup,glabel,gname,multigroup] = mgrp2idx(group,rows,sep); 
%MGRP2IDX Convert multiple grouping variables to index vector 
%   [OGROUP,GLABEL,GNAME,MULTIGROUP] = MGRP2IDX(GROUP,ROWS) takes 
%   the inputs GROUP, ROWS, and SEP.  GROUP is a grouping variable (numeric 
%   vector, string matrix, or cell array of strings) or a cell array 
%   of grouping variables.  ROWS is the number of observations. 
%   SEP is a separator for the grouping variable values. 
% 
%   The output OGROUP is a vector of group indices.  GLABEL is a cell 
%   array of group labels, each label consisting of the values of the 
%   various grouping variables separated by the characters in SEP. 
%   GNAME is a cell array containing one column per grouping variable 
%   and one row for each distinct combination of grouping variable 
%   values.  MULTIGROUP is 1 if there are multiple grouping variables 
%   or 0 if there are not. 
 
%   Tom Lane, 12-17-99 
%   Copyright 1993-2002 The MathWorks, Inc.  
%   $Revision: 1.4 $  $Date: 2002/02/04 19:25:44 $ 
 
multigroup = (iscell(group) & size(group,1)==1); 
if (~multigroup) 
   [ogroup,gname] = grp2idx(group); 
   glabel = gname; 
else 
   % Group according to each distinct combination of grouping variables 
   ngrps = size(group,2); 
   grpmat = zeros(rows,ngrps); 
   namemat = cell(1,ngrps); 
    
   % Get integer codes and names for each grouping variable 
   for j=1:ngrps 
      [g,gn] = grp2idx(group{1,j}); 
      grpmat(:,j) = g; 
      namemat{1,j} = gn; 
   end 
    
   % Find all unique combinations 
   [urows,ui,uj] = unique(grpmat,'rows'); 
    
   % Create a cell array, one col for each grouping variable value 
   % and one row for each observation 
   ogroup = uj; 
   gname = cell(size(urows)); 
   for j=1:ngrps 
      gn = namemat{1,j}; 
      gname(:,j) = gn(urows(:,j)); 
   end 
    
   % Create another cell array of multi-line texts to use as labels 
   glabel = cell(size(gname,1),1); 
   if (nargin > 2) 
      nl = sprintf(sep); 
   else 
      nl = sprintf('\n'); 
   end 
   fmt = sprintf('%%s%s',nl); 
   lnl = length(fmt)-3;        % one less than the length of nl 
   for j=1:length(glabel) 
      gn = sprintf(fmt, gname{j,:}); 
      gn(end-lnl:end) = []; 
      glabel{j,1} = gn; 
   end 
end 
end

function circleCO = circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
circleCO = [x+xp;y+yp];
end
