load 'g3.mat'
route = ShortestPathDijkstra (g3(:, 1: 2), g3(:, 3), 1, 256)