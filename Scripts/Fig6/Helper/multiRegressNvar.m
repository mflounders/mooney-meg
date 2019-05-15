function rsquared = multiRegressNvar(DV, IVN)
%Multi regression (spearman)

DV = tiedrank(DV,0);
numVars = size(IVN,1);
IVN_ranked = zeros(size(IVN));
for i_var = 1:numVars
    IVN_ranked(i_var,:) = tiedrank(IVN(i_var,:),0);
end

tableData = table;
tableData.a = DV;

analysisString = 'a~';
for i_var = 1:numVars-1
    tableData.(['b' num2str(i_var)]) = IVN_ranked(i_var,:)';
    analysisString = [analysisString 'b' num2str(i_var) '+'];
end
tableData.(['b' num2str(numVars)]) = IVN_ranked(numVars,:)';
analysisString = [analysisString 'b' num2str(numVars)];
linearModel = fitlm(tableData,analysisString);                                                      
rsquared = linearModel.Rsquared.Ordinary;