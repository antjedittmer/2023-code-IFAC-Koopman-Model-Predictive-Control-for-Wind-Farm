function writeWRCModelTable(VAFall, resCL, featureSets, modelNames, maxEig, fileOut)
% writeWRCModelTable writes the comparison of the input-only Koopman power
% models of the WRC-AIC chapter as LaTeX table (booktabs): per model (row)
% the VAF of Ur1, Ur2 and the farm power on the diagnostic data sets, the
% largest eigenvalue magnitude of A and TE / AA of the three closed loops
%
% Inputs
% - VAFall: table of plotPowerStatesDiag with column Model (feature set)
% - resCL: closed loop results (columns Model, Case, TE_kW, AA from 500 s,
%   optional TE240_kW, AA240 from sample 240 on, in brackets)
% - featureSets: feature sets of the rows (cell, order of the rows)
% - modelNames: names of the rows (cell)
% - maxEig: max|eig(A)| per row
% - fileOut: file name of the table (.tex)

sets = unique(VAFall.Set(VAFall.Model == featureSets{1}), 'stable');
types = VAFall.Type(VAFall.Model == featureSets{1});
cases = unique(resCL.Case, 'stable');

fid = fopen(fileOut, 'w');
fprintf(fid, '%% Created by writeWRCModelTable (mainGeneratePlotThesis) on %s\n', char(datetime('today','Format','yyyy-MM-dd')));
fprintf(fid, '%s\n', '\begin{table}[h]');
fprintf(fid, '%s\n', '    \centering');
fprintf(fid, '%s\n', ['    \caption[Input-only Koopman power models of the WRC-AIC]{\label{tab:WRCModelComparison} ', ...
    'Input-only Koopman power models of the WRC-AIC: VAF of the effective wind speeds $V_{r1}$, $V_{r2}$ and ', ...
    'the wind farm power $P_\text{WF}$ (simulated from the measured initial state, VAF values below 0 are set to 0), ', ...
    'largest eigenvalue magnitude of $A$ and closed-loop tracking error TE and actuator activity AA ', ...
    'from 500\,s on (as in Figures~\ref{fig:CLbaseline} to~\ref{fig:CLMPC2}) and, in brackets, from sample 240 on.}']);
fprintf(fid, '%s\n', '    \resizebox{\textwidth}{!}{%');
fprintf(fid, '    \\begin{tabular}{@{}l*{%d}{c}@{}}\n', numel(sets) + 1 + numel(cases));
fprintf(fid, '%s\n', '        \toprule');
fprintf(fid, '         & \\multicolumn{%d}{c}{VAF $V_{r1}$ / $V_{r2}$ / $P_\\text{WF}$ (\\%%)} & & \\multicolumn{%d}{c}{TE (kW) / AA $\\cdot 10^{-3}$} \\\\\n', ...
    numel(sets), numel(cases));
fprintf(fid, '        \\cmidrule(lr){2-%d}\\cmidrule(lr){%d-%d}\n', numel(sets) + 1, numel(sets) + 3, numel(sets) + 2 + numel(cases));
head = 'Model';
for idx = 1:numel(sets)
    head = sprintf('%s & %s \\texttt{%s}', head, types(idx), strrep(sets(idx),'_','\_'));
end
head = [head, ' & $\max|\lambda(A)|$'];
for idx = 1:numel(cases)
    head = sprintf('%s & %s', head, cases(idx));
end
fprintf(fid, '        %s \\\\\n', head);
fprintf(fid, '%s\n', '        \midrule');
for idxR = 1:numel(featureSets)
    row = modelNames{idxR};
    for idx = 1:numel(sets)
        v = VAFall(VAFall.Model == featureSets{idxR} & VAFall.Set == sets(idx), :);
        row = sprintf('%s & %.1f / %.1f / %.1f', row, v.VAF_Ur1, v.VAF_Ur2, v.VAF_P);
    end
    row = sprintf('%s & %.4f', row, maxEig(idxR));
    for idx = 1:numel(cases)
        r = resCL(resCL.Model == featureSets{idxR} & resCL.Case == cases(idx), :);
        row = sprintf('%s & %.1f / %.1f', row, r.TE_kW, r.AA*1e3);
        if ismember('TE240_kW', r.Properties.VariableNames) && ~isnan(r.TE240_kW)
            row = sprintf('%s (%.1f / %.1f)', row, r.TE240_kW, r.AA240*1e3);
        end
    end
    fprintf(fid, '        %s \\\\\n', row);
end
fprintf(fid, '%s\n', '        \bottomrule');
fprintf(fid, '%s\n', '    \end{tabular}}');
fprintf(fid, '%s\n', '\end{table}');
fclose(fid);
fprintf('Table written: %s\n', fileOut);
end
