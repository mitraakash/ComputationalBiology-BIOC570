function [ output_args ] = CorrelationMatrix(matrix)

n=size(matrix,1);         %Determine number of rows with "size" and"1"
standard=zscore(matrix);  %Standardize every column in matrix with "zscore"
x=2:8;                  %Set first range of columns for correlation
y=9:12;               %Set second range of columns for correlation
%%Use matrix multiplication to multiply column 1 by columns 926:1015, then
%%column 2 by columns 926:1015, and so on until all multiplications are
%%complete.  Matrix multiplication sums the values in every column
%%automatically, and this time it leaves one row vector of length 90 for
%%each of the 925 separate loops.  This equation completes the equations
%%for finding the correlation coefficient:
coefficients=(standard(:,x)'*standard(:,y))/(n-1);
%Save correlation variable to MAT file:
save('MatrixCoefficients', 'coefficients')
%Save correlation variable "correlationResults" to XLSX file:
xlswrite('MatrixCoefficients.xlsx', coefficients)
end