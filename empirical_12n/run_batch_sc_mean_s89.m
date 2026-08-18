% Batch: SC binary + graded RegLin-Mean 12n for server 89 subjects (12,14,16,18).
% Run AFTER reglin_mean iter2 for S16 and S18 complete.
for s = [12, 14, 16, 18]
    fprintf('\n=== SC-binary S%d ===\n', s);
    run_sc_binary_reglin_mean12n(s);
    fprintf('\n=== SC-graded S%d ===\n', s);
    run_sc_graded_reglin_mean12n(s);
end
fprintf('\nAll SC mean DCMs done for S89 subjects.\n');
