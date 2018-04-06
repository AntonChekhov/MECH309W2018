
loglog(nArray, ...
    abs(convergenceArr-convergenceArr(length(convergenceArr))))
% plot(nArray, ...
%     abs(convergenceArr-convergenceArr(length(convergenceArr))))
hold on 
loglog(nArray) 

nTest = log(nArray);
convTest = log(abs(convergenceArr-convergenceArr(length(convergenceArr))));
compSeries = arrayfun(@(n) 1/n^1.7, nArray); compSeries = log(compSeries)-5;

figure
plot(nTest, convTest);
hold on
plot(nTest, compSeries);
legend('Convergence Error', '1/n^{1.7}')
