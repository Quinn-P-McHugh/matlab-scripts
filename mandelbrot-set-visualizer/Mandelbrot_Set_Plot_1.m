%% Mandelbrot Set Plot Generator
% Plots the Mandelbrot set in the range (x? -2 : 1) (y? -1 : 1).  The plot is 601 x 401 pixels in size.

% Quinn McHugh
% Mechanical Engineering Lab - Section 1
% Professor Constans
% 3/29/2017

%% Initialize Variables
clear variables; close all; clc

maxIterations = 200;    % Sets the max number of iterations for each pixel

xPixels = 601;   % Sets the number of pixels along the x-axis
yPixels = 401;   % Sets the number of pixels along the y-axis

xLowerBound = -2;   % Sets the lower x bound of the grid
xHigherBound = 1;   % Sets the higher x bound of the grid
yLowerBound = -1;   % Sets the lower y bound of the grid
yHigherBound = 1;   % Sets the high y bound of the grid

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
