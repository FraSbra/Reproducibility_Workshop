clear;
close all;
clc;

saveResults_fig = '../../results/fig';
saveResults_tab = '../../results/tab';
addpath(genpath("../../excel"));
addpath(genpath("../routines"));

run("data_build.m"); % Construct multiple pred quantiles for several specific
run("table1.m"); % Test in-sample accuracy 
run("table2.m"); % Test OOS accuracy 
run("figure1.m"); % Down-Up uncertainty
