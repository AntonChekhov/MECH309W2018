
loglog(nArray, ...
    abs(convergenceArr-convergenceArr(length(convergenceArr))))
% plot(nArray, ...
%     abs(convergenceArr-convergenceArr(length(convergenceArr))))


nTest = log(nArray);
convTest = log(abs(convergenceArr-convergenceArr(length(convergenceArr))));

figure
plot(nTest, convTest);
convTest(16)/nTest(16)
