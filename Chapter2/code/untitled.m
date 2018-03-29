load('data/Processed/webGeneration.mat')

bodySizesFig = figure('Units','inches','Position',[0 0 12 12]);

ax2 = subplot(2,1,2);
nSpeciesWithBS = cellfun(@(x) isfinite(x(:,12)),propertiesCell(1:3),'UniformOutput',false);
allConBS = cell2mat(cellfun(@(x,y) x(y(:,2),12),propertiesCell(1:3,1),linkListCell(1:3),'UniformOutput',false));
allResBS = cell2mat(cellfun(@(x,y) x(y(:,1),12),propertiesCell(1:3,1),linkListCell(1:3),'UniformOutput',false));
goodBS = isfinite(allConBS)&isfinite(allResBS);
allConBS = allConBS(goodBS);
allResBS = allResBS(goodBS);
linkTypes = cell2mat(cellfun(@(x,y,z) x(y(:,2)) + 3*z(y(:,2)),speciesTypeCell(1:3),linkListCell(1:3),propertiesCell(1:3,2),'UniformOutput',false));
linkTypes = linkTypes(goodBS);
ptColors = [0 0.7 0;0 0 0.7; 0 0 0;0.7 0 0];
gscatter(log10(allResBS),log10(allConBS),linkTypes,ptColors,'.',8)
r2s = (grpstats([log10(allConBS),log10(allResBS)],linkTypes,@corr)).^2;
R2s = r2s(:,2,1);
xl2 = xlim;
yl2 = ylim;

title('(b) Consumer and Resource Body Sizes by Consumer Type in Trophic Webs','FontName','CMU Serif')
ylabel('$\log_{10}($Con. Body Size$)$','Interpreter','LaTeX');
xlabel('$\log_{10}($Res. Body Size$)$','Interpreter','LaTeX');
rl = refline(1,0);
rl.LineStyle = '--';
rl.Color = [0.7 0.7 0.7];


lsl = lsline;
slopes = arrayfun(@(x) diff(x.YData)/diff(x.XData),lsl);
intercepts = arrayfun(@(x,y) y.YData(1) - x*y.XData(1) ,slopes,lsl);

leg = legend;
leg.FontName = 'CMU Serif';
leg.Interpreter = 'LaTeX';
leg.Orientation = 'horizontal';
leg.Location = 'south';
legend(sprintf('Invertebrates; $R^2=%.2f$\n$y=%.2fx %+.2f$',R2s(1),slopes(4),intercepts(4))...
            ,sprintf('Ectotherm Vert.; $R^2=%.2f$\n$y=%.2fx %+.2f$',R2s(2),slopes(3),intercepts(3))...
            ,sprintf('Endotherm; $R^2=%.2f$\n$y=%.2fx %+.2f$',R2s(3),slopes(2),intercepts(2))...
            ,sprintf('Parasites; $R^2=%.2f$\n$y=%.2fx %+.2f$',R2s(4),slopes(1),intercepts(1))...
            );

ax1 = subplot(2,1,1);
gscatter(log10(cell2mat(OGResBSCell)),log10(cell2mat(OGConBSCell))...
         ,cell2mat(OGLinkTypeCell)...
         ,ptColors,'.',8)
r2s = (grpstats([log10(cell2mat(OGResBSCell)),log10(cell2mat(OGConBSCell))]...
         ,cell2mat(OGLinkTypeCell),@corr)).^2;
R2s = r2s(:,2,1);
leg = legend('Invertebrates'...
            ,'Ectotherm Vert.'...
            ,'Endotherm'...r
            ,'Parasites');
leg.FontName = 'CMU Serif';

title('(a) Consumer and Resource Body Sizes by Consumer Type in Unaggregated Webs','FontName','CMU Serif')
ylabel('$\log_{10}($Con. Body Size$)$','Interpreter','LaTeX');
xlabel('$\log_{10}($Res. Body Size$)$','Interpreter','LaTeX');
rl = refline(1,0);
rl.LineStyle = '--';
rl.Color = [0.7 0.7 0.7];

lsl = lsline;
slopes = arrayfun(@(x) diff(x.YData)/diff(x.XData),lsl);
intercepts = arrayfun(@(x,y) y.YData(1) - x*y.XData(1), slopes,lsl);

leg = legend;
leg.FontName = 'CMU Serif';
leg.Interpreter = 'LaTeX';
leg.Orientation = 'horizontal';
leg.Location = 'south';
legend(sprintf('Invertebrates; $R^2=%.2f$\n$y=%.2f x %+.2f$',R2s(1),slopes(4),intercepts(4))...
            ,sprintf('Ectotherm Vert.; $R^2=%.2f$\n$y=%.2f x %+.2f$',R2s(2),slopes(3),intercepts(3))...
            ,sprintf('Endotherm; $R^2=%.2f$\n$y=%.2f x %+.2f$',R2s(3),slopes(2),intercepts(2))...
            ,sprintf('Parasites; $R^2=%.2f$\n$y=%.2f x %+.2f$',R2s(4),slopes(1),intercepts(1))...
            );
arrayfun(@(x) set(x,'FontName','CMU Serif'),bodySizesFig.Children)
xl1 = xlim;
yl1 = ylim;

xl = [min(xl1(1),xl2(1)) max(xl1(2),xl2(2))];
yl = [min(yl1(1),yl2(1)) max(yl1(2),yl2(2))];
xlim(ax1,xl);
ylim(ax1,yl);

xlim(ax2,xl);
ylim(ax2,yl);


print(sprintf('../figures/bodySizes.png'),'-dpng','-r0')
