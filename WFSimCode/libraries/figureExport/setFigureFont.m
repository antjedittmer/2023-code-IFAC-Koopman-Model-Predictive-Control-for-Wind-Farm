function setFigureFont(h, fs)
%SETFIGUREFONT One font size for all text of a figure, bold titles
%   setFigureFont(h) sets the font size of all axes, labels, titles,
%   legends and text objects of figure h to 12 pt (as the thesis figures)
%   and the titles to bold. Call it after all text is created.
%   setFigureFont(h, fs) uses the font size fs.

if nargin < 2 || isempty(fs)
    fs = 12;
end
set(findall(h, '-property', 'FontSize'), 'FontSize', fs);
axs = findall(h, 'Type', 'axes');
set(axs, 'TitleFontSizeMultiplier', 1, 'LabelFontSizeMultiplier', 1);
for idx = 1:numel(axs)
    axs(idx).Title.FontWeight = 'bold';
    axs(idx).Title.FontSize = fs;
end
end
