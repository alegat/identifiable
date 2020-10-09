% Tests 6 small networks taken from the literature.

% In main_identifiable(name, graph, excited, measured, nsamples),
% nsamples = 1 should be enough since problematic cases lie in a zero-
% measure set, but we take nsamples = 5 to be sure to avoid numerical
% issues.

function test_identifiable
close all;

% Fig 1 HGB18
graph = [0 0 0;
         1 0 1;
         1 0 0];
excited = [1 1 1]';
measured = [0 1 1]'; % identifiable
name = 'Fig 1 HGB18';
main_identifiable(name, graph, excited, measured, 0);

% Fig 3 HGB18
graph = [0 1 1;
         1 0 1;
         0 1 0];
excited = [1 1 1]';
%measured = [1 1 0]'; % identifiable
measured = [0 0 1]'; % Could not identify 4 edges
name = 'Fig 3 HGB18';
main_identifiable(name, graph, excited, measured, 5);

% Ex 1 HGB18
graph = [0 0 0 0 0;
         1 0 0 0 0;
         1 0 0 0 0;
         0 1 1 0 0;
         0 1 1 0 0];
excited = [1 1 1 1 1]';
%measured = [0 0 0 1 1]'; % identifiable
measured = [0 0 1 0 1]'; % Could not identify 2 edges
name = 'Ex. 1 HGB18';
main_identifiable(name, graph, excited, measured, 5);

% Fig 6/8 CSV19
%        1 2 3 4 5 6 7 8 9 1011
graph = [0 1 0 0 0 0 0 0 0 0 0; % 1
         0 0 1 0 0 0 0 0 0 0 0; % 2
         0 1 0 0 0 0 0 0 0 0 0; % 3
         0 0 0 0 0 0 0 0 0 0 0; % 4
         1 0 0 0 0 1 0 0 1 0 0; % 5
         0 1 0 0 0 0 1 0 1 0 0; % 6
         0 0 1 1 0 0 0 0 0 1 1; % 7
         0 0 0 1 0 0 1 0 0 0 1; % 8
         0 0 0 0 0 0 1 0 0 1 0; % 9
         0 0 0 0 0 1 0 0 0 0 0; % 10
         0 0 0 0 0 0 0 0 0 0 0];% 11
%measured = ones(11,1);
%excited = [0 1 1 1 0 0 0 0 0 0 1]'; % Could not identify 6 edges
measured = [1 1 1 0 1 1 1 1 1 0 0]';
excited =  [0 1 0 1 0 0 0 0 0 1 1]'; % identifiable
name = 'Fig 8 CSV19';
main_identifiable(name, graph, excited, measured, 5);

% Fig 9 CSV19
%        1 2 3 4 5 6 7 8 9
graph = [0 0 0 0 0 0 0 0 0; % 1
         1 0 0 1 0 0 0 0 0; % 2
         0 0 0 0 1 0 0 0 0; % 3
         0 0 0 0 0 0 0 0 0; % 4
         0 1 0 1 0 0 0 1 0; % 5
         0 0 0 0 1 0 0 1 0; % 6
         0 0 0 0 1 0 0 0 0; % 7
         0 0 0 0 0 0 1 0 0; % 8
         0 0 0 0 0 0 0 1 0]; % 9
excited = ones(9,1);
%measured = [0 1 1 0 0 1 0 0 1]'; % identifiable
measured = [0 1 1 0 0 1 0 1 0]'; % Could not identify 1 edge
%measured = [0 1 1 1 0 0 0 1 0]'; % Could not identify 3 edges
name = 'Fig 9 CSV19';
main_identifiable(name, graph, excited, measured, 5);

end
