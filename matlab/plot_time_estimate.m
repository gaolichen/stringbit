function f = plot_time_estimate(bits, statenumber, pointNumber)
    tic;
    states = lowest_energies(1, bits, bits, statenumber);
    elapsedTime = toc;
    f = elapsedTime * pointNumber;
end