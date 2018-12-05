function graph_output = chainRule(graph)
%ChainRule is a function which accumulated derivatives in reverse order, 
% from end node to multiple start nodes 
% "graph" is the top-node of the calculation tree. Call this top-node
% "current" for simplicity:
current = graph;

% Set derivative to 1: (df/df = 1)
current.derivative = 1;


% Push derivative to parents, and grandparents etc.. :
push_derivative(current);

% Return modified graph:
graph_output = graph;
end

function push_derivative(node)
    if isempty(node.parents)
        return
    end
    % Push derivative from current node to the nodes parents, according to
    % corresponding derivativeOp-function:
    node.derivativeOp(node.derivative, node.parents);
    % Since there exist at least one parent, this is the left one. Call the push_derivative function
    % on the current nodes left-parent.
    push_derivative(node.parents(1))
    % If node has "right-parent", then call the push_derivative function
    % on the current nodes right-parent.
    if size(node.parents,2)==2
        push_derivative(node.parents(2))
    end
    % We have now pused all derivatives down to all start nodes. 
end



