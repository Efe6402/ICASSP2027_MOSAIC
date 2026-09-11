function generate_coil20_paper_figures(fit_file,data_file,output_dir)
%GENERATE_COIL20_PAPER_FIGURES Save the four six-view COIL-20 paper panels.

if ~isfile(fit_file), error('Fit file not found: %s',fit_file); end
if ~isfile(data_file), error('Data file not found: %s',data_file); end
if ~isfolder(output_dir), mkdir(output_dir); end

S = load(fit_file);
if isfield(S,'compactFit'), fit=S.compactFit;
elseif isfield(S,'primaryFit'), fit=S.primaryFit;
elseif isfield(S,'fit'), fit=S.fit;
else, error('The fit file does not contain a fitted MOSAIC model.');
end

D=load(data_file,'data');
if ~isfield(D,'data'), error('The data file does not contain data.'); end
data=normalize_data(D.data);

if isfield(S,'analysis') && isfield(S.analysis,'within_view_adjacencies')
    within=S.analysis.within_view_adjacencies;
elseif isfield(S,'primaryAnalysis') && ...
        isfield(S.primaryAnalysis,'within_view_adjacencies')
    within=S.primaryAnalysis.within_view_adjacencies;
else
    within=within_view_adjacencies(fit,data.num_views);
end

plot_modal_profiles(fit.B,data.view_angles_deg, ...
    fullfile(output_dir,'01_modal_view_profiles.png'));
plot_within_graph(data,within(:,:,1),1, ...
    fullfile(output_dir,'02_within_view_graph.png'));

edges=select_diverse_cross_edges(fit.A_supra,data.num_objects, ...
    data.num_views,12,1,3);
plot_cross_pairs(data,edges, ...
    fullfile(output_dir,'03_crossnode_crossview_pairs.png'));
plot_view_graph(fit.view_graph,data.view_angles_deg, ...
    fullfile(output_dir,'04_induced_view_graph.png'));
end

function data=normalize_data(data)
data.num_views=numel(data.X);
data.num_objects=size(data.X{1},1);
data.object_ids=double(data.object_ids(:));
data.view_angles_deg=double(data.view_angles_deg(:)');
if isfield(data,'processed_image_size')
    data.processed_image_size=double(data.processed_image_size(:)');
else
    side=round(sqrt(size(data.X{1},2)));
    data.processed_image_size=[side side];
end
if data.num_views~=6, error('Expected six views.'); end
end

function within=within_view_adjacencies(fit,K)
n=size(fit.A_state,1); r=size(fit.B,2); within=zeros(n,n,K);
for k=1:K
    for m=1:r
        within(:,:,k)=within(:,:,k)+fit.B(k,m)^2*fit.A_state(:,:,m);
    end
end
end

function plot_modal_profiles(B,angles,file_path)
r=size(B,2);
f=figure('Visible','off','Color','w','Position',[30 30 1750 520]);
tl=tiledlayout(f,1,r,'TileSpacing','compact','Padding','compact');
for m=1:r
    ax=nexttile(tl,m);
    plot(ax,[angles 360],[B(:,m);B(1,m)],'-o','LineWidth',1.7, ...
        'MarkerSize',6,'MarkerFaceColor','w');
    grid(ax,'on'); box(ax,'on'); xlim(ax,[0 360]);
    ylim(ax,[0 max(0.6,1.08*max(B(:,m)))]);
    xticks(ax,0:60:360);
    xlabel(ax,'angle (deg)'); ylabel(ax,'b_m weight');
    title(ax,sprintf('Mode %d view profile',m),'FontWeight','bold');
end
export_png(f,file_path,300); close(f);
end

function plot_within_graph(data,A,view_index,file_path)
n=data.num_objects;
[ei,ej,ew]=top_undirected_edges(A,15);
f=figure('Visible','off','Color','w','Position',[30 30 1550 1450]);
ax=axes('Parent',f,'Position',[0.035 0.035 0.93 0.93]);
hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
xlim(ax,[-1.52 1.52]); ylim(ax,[-1.52 1.52]); set(ax,'YDir','normal');
colormap(ax,gray(256));
theta=pi/2-2*pi*(0:n-1)/n;
centers=1.15*[cos(theta(:)),sin(theta(:))];
for q=1:numel(ew)
    p1=centers(ei(q),:); p2=centers(ej(q),:);
    [s,t]=segment_boundaries(p1,p2,.095,.095);
    plot(ax,[s(1) t(1)],[s(2) t(2)],'-','Color',[.08 .54 .70], ...
        'LineWidth',.6+4.7*ew(q)/max(ew));
end
for i=1:n
    image_matrix=reshape(data.X{view_index}(i,:),data.processed_image_size);
    draw_image(ax,image_matrix,centers(i,1),centers(i,2),.19,.19);
    radius=1.305;
    text(ax,radius*cos(theta(i)),radius*sin(theta(i)), ...
        sprintf('%d',data.object_ids(i)),'HorizontalAlignment','center', ...
        'VerticalAlignment','middle','FontSize',8,'FontWeight','bold');
end
export_png(f,file_path,300); close(f);
end

function edges=select_diverse_cross_edges(A,n,K,maximum_edges,max_object,max_view)
template=struct('i',NaN,'j',NaN,'k',NaN,'ell',NaN,'weight',NaN);
records=repmat(template,0,1);
for k=1:K
    Ik=(k-1)*n+(1:n);
    for ell=k+1:K
        Il=(ell-1)*n+(1:n); block=A(Ik,Il);
        for i=1:n
            for j=1:n
                if i~=j && isfinite(block(i,j)) && block(i,j)>0
                    z=template; z.i=i; z.j=j; z.k=k; z.ell=ell;
                    z.weight=double(block(i,j)); records(end+1,1)=z; %#ok<AGROW>
                end
            end
        end
    end
end
[~,order]=sort([records.weight],'descend'); records=records(order);
object_count=zeros(n); view_count=zeros(K); edges=records([]);
for q=1:numel(records)
    z=records(q); ia=min(z.i,z.j); ib=max(z.i,z.j);
    ka=min(z.k,z.ell); kb=max(z.k,z.ell);
    if object_count(ia,ib)>=max_object || view_count(ka,kb)>=max_view
        continue;
    end
    edges(end+1,1)=z; %#ok<AGROW>
    object_count(ia,ib)=object_count(ia,ib)+1;
    view_count(ka,kb)=view_count(ka,kb)+1;
    if numel(edges)>=maximum_edges, break; end
end
end

function plot_cross_pairs(data,edges,file_path)
if isempty(edges), error('No positive cross-node/cross-view edges found.'); end
number_rows=numel(edges);
f=figure('Visible','off','Color','w', ...
    'Position',[30 30 1800 max(1500,175*number_rows+150)]);
ax=axes('Parent',f,'Position',[.025 .025 .95 .95]);
hold(ax,'on'); axis(ax,'equal'); xlim(ax,[0 1]); ylim(ax,[0 1]);
axis(ax,'off'); set(ax,'YDir','normal'); colormap(ax,gray(256));
row_y=linspace(.955,.045,number_rows); row_spacing=abs(row_y(1)-row_y(2));
image_side=.82*row_spacing; image_half=image_side/2;
horizontal_gap=.14;
left_x=.5-.5*(horizontal_gap+image_side);
right_x=.5+.5*(horizontal_gap+image_side);
left_inner=left_x+image_half; right_inner=right_x-image_half;
left_outer=left_x-image_half; right_outer=right_x+image_half;
maximum_weight=max([edges.weight]);
for q=1:number_rows
    z=edges(q); y=row_y(q);
    plot(ax,[left_inner+.006 right_inner-.006],[y y],'-', ...
        'Color',[.08 .54 .70],'LineWidth',1+6*z.weight/maximum_weight);
    left_image=reshape(data.X{z.k}(z.i,:),data.processed_image_size);
    right_image=reshape(data.X{z.ell}(z.j,:),data.processed_image_size);
    draw_image(ax,left_image,left_x,y,image_side,image_side);
    draw_image(ax,right_image,right_x,y,image_side,image_side);
    text(ax,left_outer-.015,y,sprintf('obj %d, %g deg', ...
        data.object_ids(z.i),data.view_angles_deg(z.k)), ...
        'HorizontalAlignment','right','VerticalAlignment','middle', ...
        'FontSize',11,'FontWeight','bold');
    text(ax,right_outer+.015,y,sprintf('obj %d, %g deg', ...
        data.object_ids(z.j),data.view_angles_deg(z.ell)), ...
        'HorizontalAlignment','left','VerticalAlignment','middle', ...
        'FontSize',11,'FontWeight','bold');
end
export_png(f,file_path,300); close(f);
end

function plot_view_graph(G,angles,file_path)
K=numel(angles); [ei,ej,ew]=top_undirected_edges(G,3);
f=figure('Visible','off','Color','w','Position',[30 30 800 800]);
ax=axes('Parent',f,'Position',[.05 .05 .9 .9]); hold(ax,'on');
axis(ax,'equal'); axis(ax,'off'); xlim(ax,[-1.35 1.35]); ylim(ax,[-1.35 1.35]);
theta=pi/2-2*pi*(0:K-1)/K; xy=[cos(theta(:)),sin(theta(:))];
for q=1:numel(ew)
    plot(ax,xy([ei(q) ej(q)],1),xy([ei(q) ej(q)],2),'-', ...
        'Color',[.25 .25 .25],'LineWidth',5);
end
colors=[0 .447 .741; .929 .286 .027; .929 .694 .125; ...
        0 .447 .741; .929 .286 .027; .929 .694 .125];
for k=1:K
    scatter(ax,xy(k,1),xy(k,2),260,colors(k,:),'filled', ...
        'MarkerEdgeColor','k','LineWidth',1.4);
    label_xy=1.18*xy(k,:);
    text(ax,label_xy(1),label_xy(2),sprintf('%g^\\circ',angles(k)), ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'Color',colors(k,:),'FontSize',14,'FontWeight','bold', ...
        'Interpreter','tex');
end
export_png(f,file_path,300); close(f);
end

function [i,j,w]=top_undirected_edges(A,maximum_edges)
mask=triu(true(size(A)),1); values=A(mask); [rows,cols]=find(mask);
keep=isfinite(values)&values>0; values=values(keep); rows=rows(keep); cols=cols(keep);
[w,order]=sort(values,'descend'); take=min(maximum_edges,numel(order));
order=order(1:take); i=rows(order); j=cols(order); w=w(1:take);
end

function draw_image(ax,image_matrix,cx,cy,width,height)
x=[cx-width/2 cx+width/2]; y=[cy+height/2 cy-height/2];
image(ax,'XData',x,'YData',y,'CData',double(image_matrix),'CDataMapping','scaled');
rectangle(ax,'Position',[cx-width/2 cy-height/2 width height], ...
    'EdgeColor',[.12 .12 .12],'LineWidth',.9);
end

function [start_point,end_point]=segment_boundaries(c1,c2,half_width,half_height)
direction=c2-c1; distance=norm(direction);
if distance<=eps, start_point=c1; end_point=c2; return; end
u=direction/distance;
offset=min([half_width/max(abs(u(1)),eps),half_height/max(abs(u(2)),eps)]);
start_point=c1+offset*u; end_point=c2-offset*u;
end

function export_png(f,file_path,resolution)
try
    exportgraphics(f,file_path,'Resolution',resolution);
catch
    print(f,file_path,'-dpng',sprintf('-r%d',resolution));
end
end
