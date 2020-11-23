function plot_results(title_str, stages, figure_num, eta, d_eta, I, d_I)
figure(figure_num);
errorbar(eta, d_eta, '-o', 'MarkerFaceColor', 'b');
ylabel('\eta');
ylim([0 inf]);

yyaxis right;
errorbar(I, d_I, '-s', 'MarkerFaceColor', 'r');
ylabel('I_{max}');
ylim([0 inf]);

title(title_str);
xlabel('Treatment stage');
xticks(1:length(stages));
xticklabels(stages);
xlim([0.75 length(stages)+0.25]);


set(gcf, 'Position', [500, 500, 400, 300])
end