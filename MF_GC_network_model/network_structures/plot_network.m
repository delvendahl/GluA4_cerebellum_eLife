% Visualize MF-GC network, i.e., Fig. 1a

% Parameters to plot
N_syn = 4; % number inputs

load(strcat('GCLconnectivity_',int2str(N_syn),'.mat'))

[N_mf,N_grc] = size(conn_mat);

dist = @(x,y) sqrt(sum((x-y).^2));
%%
plot_1 = figure(1, 'position',[100,100,500,500]);
hold on;
% Plot connections and get dendritic lengths
dend = [];
for k1 = 1:N_mf
    for k2 = 1:N_grc
        if conn_mat(k1,k2)==1
            h = plot3([glom_pos(k1,1),grc_pos(k2,1)],[glom_pos(k1,2),grc_pos(k2,2)],[glom_pos(k1,3),grc_pos(k2,3)],'Color',[0.5,0.5,0.5],'LineWidth',3);
            dend = [dend, dist(glom_pos(k1,:),grc_pos(k2,:))];
        end
    end
end


[x,y,z] = sphere; r = 2;

% Plot MF glomeruli
for k1 = 1:N_mf
    h=surfl(r*x+glom_pos(k1,1),r*y+glom_pos(k1,2),r*z+glom_pos(k1,3));
    set(h,'FaceColor',[1, 0, 0]);
    set(h,'EdgeColor','none');
    set(h,'FaceLighting', 'Gouraud');
    set(h,"DiffuseStrength", 0.1*(1+rand/3), "AmbientStrength", .5*(1+rand/5), "SpecularStrength", 0.3*(1+rand/3));
end
light('Position',[-100 -50 0],'Style','infinite');
    
% Plot GCs
for k1 = 1:N_grc
    h=surfl(r*x+grc_pos(k1,1),r*y+grc_pos(k1,2),r*z+grc_pos(k1,3));
    set(h,'FaceColor',[0, 0, 1]);
    set(h,'EdgeColor','none');
    set(h,'FaceLighting', 'Gouraud');
    set(h,"DiffuseStrength", 0.1*(1+rand/3), "AmbientStrength", .5*(1+rand/5), "SpecularStrength", 0.3*(1+rand/3));
end

light('Position',[-100 -50 0],'Style','infinite');
    
print(plot_1, "fig1.pdf")
print('fig1.png', '-dpng', '-S1200,1200');

%% Plot distribution of dendritic lengths
bins = 0:2.5:40; h = histc(dend,bins);
plot_2 = figure;
b = bar(bins,h); xlim([0,42])
hold on, plot(median(dend),950,'v','Color',[.35,0,.5])
set(gca,'FontSize',20)
xlabel('Dendritic length (\mum)')
ylabel('Number')
set(b,'EdgeColor','w','FaceColor',[.35,0,.5])

print(plot_2, "fig2.pdf")