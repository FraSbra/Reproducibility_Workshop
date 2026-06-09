clear;
close all;
clc;

saveResults_fig = '../../results/fig';
saveResults_tab = '../../results/tab';
addpath(genpath("../../excel")); % Folder containing excel files
addpath(genpath("../routines")); % Folder contains all functions


run("figure1.m"); % IRF of eonia for US MP shock
run("figure2.m"); % IRF of risk for US MP shock 
run("figure3.m"); % Hold fix EU MP response
run("figure4.m"); % Hold fix financial market response
run("figure5.m"); % IRF of risk for EA MP shock
