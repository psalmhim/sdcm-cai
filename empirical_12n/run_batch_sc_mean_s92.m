% Batch: SC binary + graded RegLin-Mean 12n for server 92 subjects (13,15,17).
% Run AFTER reglin_mean iter2 for S13, S15, S17 complete.
for s = [13, 15, 17]
    fprintf('\n=== SC-binary S%d ===\n', s);
    run_sc_binary_reglin_mean12n(s);
    fprintf('\n=== SC-graded S%d ===\n', s);
    run_sc_graded_reglin_mean12n(s);
end
fprintf('\nAll SC mean DCMs done for S92 subjects.\n');
