h = openfig('D15_ALC_all_thetas.fig','invisible');
ax = findobj(h,'Type','axes');
ln = findobj(ax,'Type','line');

x = get(ln,'XData');
y = get(ln,'YData');

thetas = linspace(0,90,200);   % degrees

n = 10;
idx = 1:n:numel(ln);
cmap = lines(numel(idx));

figure; hold on
for i = 1:numel(idx)
    k = idx(i);
    plot(x{k}, y{k}, 'Color', cmap(i,:));
end
hold off

legendStrings = arrayfun(@(k) sprintf('\\theta = %.1f^\\circ', thetas(k)), ...
                          idx, 'UniformOutput', false);
legend(legendStrings, 'Location', 'best');
