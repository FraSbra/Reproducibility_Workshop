function [dx] = MatrixDiff(xx)

dx = xx - mat_single_leadlag(xx,-1); 
