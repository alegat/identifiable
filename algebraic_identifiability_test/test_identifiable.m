function test_identifiable
close all;

graph = [0 0 0;      % the adjacency matrix
         1 0 1;
         1 0 0];
unknowns = [0 0 0;   % the unknown edges
            1 0 0;
            1 0 0];
excited = [1 0 0]';  % the excited nodes
measured = [0 1 1]'; % the measured nodes
name = 'Example A';
nsamples = 100;
decoupled = 1;	     % Take decoupled = 1 for generic decoupled identifiability.
                     % Otherwise, it is generic local identifiability.
main_identifiable(name, graph, unknowns, excited, measured, nsamples, decoupled);

% Fig 1 HGB18
graph = [0 0 0;
         1 0 1;
         1 0 0];
unknowns = graph;
excited = [1 1 1]';
measured = [0 1 1]';
name = 'Fig 1 HGB18';
main_identifiable(name, graph, unknowns, excited, measured, 0, 0);

% Fig 3 HGB18
graph = [0 1 1;
         1 0 1;
         0 1 0];
unknowns = graph;
excited = [1 1 1]';
measured = [1 1 0]';
name = 'Fig 3 HGB18';
main_identifiable(name, graph, unknowns, excited, measured, 5, 0);

% Ex 1 HGB18
graph = [0 0 0 0 0;
         1 0 0 0 0;
         1 0 0 0 0;
         0 1 1 0 0;
         0 1 1 0 0];
unknowns = graph;
excited = [1 1 1 1 1]';
measured = [0 0 0 1 1]';
name = 'Ex 1 HGB18';
main_identifiable(name, graph, unknowns, excited, measured, 5, 0);

% Nodes: 1 2 3
graph = [0 0 1 ;
         1 0 1 ;
         1 0 0];
unknowns = [0 0 1 ;
            1 0 0 ;
            0 0 0];
measured = [0 1 1]';
excited =  [1 0 1]'; 
name = 'Example from TAC';
decoupled = 0;
main_identifiable(name, graph, unknowns, excited, measured, 10, decoupled);

decoupled = 1;
main_identifiable(name, graph, unknowns, excited, measured, 10, decoupled);
     
% Nodes: 1 2 3 4 5
graph = [0 0 1 0 1;
         1 0 0 0 0;
         0 1 0 0 0;
         0 0 1 0 0;
         1 0 0 1 0];
unknowns = [0 0 0 0 1;
            1 0 0 0 0;
            0 0 0 0 0;
            0 0 1 0 0;
            0 0 0 1 0];
measured = [1 0 0 1 0]';
excited =  [1 0 0 0 1]'; 
name = 'Example B';
main_identifiable(name, graph, unknowns, excited, measured, 10, 0);

measured = [1 0 0 0 1]';
excited =  [1 0 0 1 0]'; 
name = 'Example B bis';
main_identifiable(name, graph, unknowns, excited, measured, 10, 0);

% Counterexample to our conjecture
graph = zeros(14,14);
graph(3,1) = 1;
graph(4,1) = 1;
graph(3,2) = 1;
graph(7,2) = 1;
graph(5,3) = 1;
graph(6,3) = 1;
graph(8,4) = 1;
graph(9,5) = 1;
graph(10,6) = 1;
graph(11,7) = 1;
graph(12,8) = 1;
graph(12,9) = 1;
graph(13,12) = 1;
graph(13,10) = 1;
graph(14,12) = 1;
graph(14,11) = 1;

unknowns = zeros(14,14);
unknowns(8,4) = 1;
unknowns(9,5) = 1;
unknowns(10,6) = 1;
unknowns(11,7) = 1;

excited = zeros(14,1); excited(1,1) = 1; excited(2,1) = 1;
measured = zeros(14,1); measured(13,1) = 1; measured(14,1) = 1;
name = 'Counterex to conjecture';
main_identifiable(name, graph, unknowns, excited, measured, 100, 0);

end