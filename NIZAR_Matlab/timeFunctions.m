function [time, solutions] = timeFunctions()
    time = zeros(30, 13);
    solutions = zeros(30, 13);

    for i = 1:30
        tic;
        [~,solutions(i,1)] = NIZAROptimizer(@sphere, bounds(50,-100,100), 25, 1400);
        time(i,1) = toc;
    end
    for i = 1:30
        tic;
        [~,solutions(i,2)] = NIZAROptimizer(@quartic, bounds(50,-1.28,1.28), 25, 1400);
        time(i,2) = toc;
    end
    for i = 1:30
        tic;
        [~,solutions(i,3)] = NIZAROptimizer(@powell_sum, bounds(50,-1,1), 25, 1400);
        time(i,3) = toc;
    end
    for i = 1:30
        tic;
        [~,solutions(i,4)] = NIZAROptimizer(@sum_squares, bounds(50,-10,10), 25, 1400);
        time(i,4) = toc;
    end
    for i = 1:30
        tic;
        [~,solutions(i,5)] = NIZAROptimizer(@schwefel_2_20, bounds(50,-100,100), 25, 1400);
        time(i,5) = toc;
    end
    for i = 1:30
        tic;
        [~,solutions(i,6)] = NIZAROptimizer(@stepint, bounds(50,-5.12,5.12), 25, 1400);
        time(i,6) = toc;
    end
    for i = 1:30
        tic;
        [~,solutions(i,7)] = NIZAROptimizer(@ridge, bounds(50,-5,5), 25, 1400);
        time(i,7) = toc;
    end
    for i = 1:30
        tic;
        [~,solutions(i,8)] = NIZAROptimizer(@neumaier_N3, bounds(15,-100,100), 25, 1400);
        time(i,8) = toc;
    end
    for i = 1:30
        tic;
        [~,solutions(i,9)] = NIZAROptimizer(@ackley_N2, bounds(2,-32,32), 25, 1400);
        time(i,9) = toc;
    end
    for i = 1:30
        tic;
        [~,solutions(i,10)] = NIZAROptimizer(@shekel_10, bounds(4,0,10), 25, 1400);
        time(i,10) = toc;
    end
    for i = 1:30
        tic;
        [~,solutions(i,11)] = NIZAROptimizer(@pressure_vessel_design, [0 100; 0 100; 10 200; 10 200], 25, 1400);
        time(i,11) = toc;
    end
    for i = 1:30
        tic;
        [~,solutions(i,12)] = NIZAROptimizer(@tension_compression_spring_design, [0.05 2; 0.25 1.3; 2 15], 25, 1400);
        time(i,12) = toc;
    end
    %for i = 1:30
    %    tic;
    %    [~,solutions(i,13)] = NIZAROptimizer(@quartic, [0.05 2; 0.25 1.3; 2 15], 25, 1400);
    %    time(i,13) = toc;
    %end
end