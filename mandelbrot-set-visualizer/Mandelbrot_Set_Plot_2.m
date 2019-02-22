%% Mandelbrot Set Plot Generator
% Plots the Mandelbrot set in the range (x? -0.05 : -0.01) (y? 0.77 : 0.81).  The plot is 501 x 501 pixels in size.

% Quinn McHugh
% Mechanical Engineering Lab - Section 1
% Professor Constans
% 3/29/2017

%% Initialize Variables
clear variables; close all; clc;

maxIterations = 200;    % Sets the max number of iterations for each pixel

xPixels = 501;   % Sets the number of pixels along the x-axis
yPixels = 501;   % Sets the number of pixels along the y-axis

xLowerBound = -0.05;   % Sets the lower x bound of the grid
xHigherBound = -0.01;   % Sets the higher x bound of the grid
yLowerBound = 0.77;   % Sets the lower y bound of the grid
yHigherBound = 0.81;   % Sets the high y bound of the grid

M = zeros(yPixels,xPixels);     % Initializes the mandelbrot set with zeros
z = zeros(yPixels,xPixels);     % Initializes the z array with zeros

%% Main Program
xRange = linspace(xLowerBound, xHigherBound, xPixels);  % Generates a linearly spaced vector to serve as the range of x values in the grid
yRange = linspace(yLowerBound, yHigherBound, yPixels);  % Generates a linearly spaced vector to serve as the range of y values in the grid

[xCoord,yCoord] = meshgrid(xRange,yRange);  % Stores a 2-D grid of coordinates based on xRange and yRange

c = xCoord + yCoord*1j;     % Sets the value of c

for n = 1:maxIterations
    z = z.*z + c;                   % Calculates the current value of z for all coordinates in the grid
    isLessThanTwo = abs(z) < 2;     % An array that stores a value of 1 (true) if the current absolute value of z is less than 2
    M = M + isLessThanTwo;          % Counts the number of iterations before each coordinate in the grid "blows up"
end

%% Plot Figure
imshow(M,jet)  % Plots the Mandelbrot set and applies a color to each pixel
axis on;
axis equal;
