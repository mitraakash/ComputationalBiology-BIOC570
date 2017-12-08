function [tableTransposed] = transposeTable(tableIn)
%this function transposes a table. 
clinical_info_arr = table2cell(tableIn);
tableTransposed = cell2table(clinical_info_arr.');
tableTransposed.Properties.RowNames = tableIn.Properties.VariableNames;
end
