% Solves the symbolic equation CTB = CT'B, and outputs the solution S

function S = solve_CTB_CTbisB(graph, excited, measured)

nb_edges = sum(graph,"all");
[L,~] = size(graph);
Id = eye(L);
B = diag(excited); [~,Bcol] = find(B); B = B(:,Bcol);
C = diag(measured);[Crow,~] = find(C); C = C(Crow,:);

syms G Gbis [L L]
G = G .* graph;
Gbis = Gbis .* graph;

eqn1 = reshape(C*inv(Id-G)*B == C*inv(Id-Gbis)*B, 1, []);
eqn2 = (reshape(nonzeros(G),1,[]) ~= zeros(1,nb_edges));
eqn3 = (reshape(nonzeros(Gbis),1,[]) ~= zeros(1,nb_edges));

eqns = [eqn1, eqn2, eqn3];
vars = nonzeros([reshape(G,1,[]) reshape(Gbis,1,[])]);
S = solve(eqns, vars, 'ReturnConditions', true);

end