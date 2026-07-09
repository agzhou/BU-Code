% assign fields from a structure 's' at some index 'ind', into a vector 'x', in order
% adapted from https://stackoverflow.com/a/29958475

function x = struct2vec(s, ind)
    
    fn = fieldnames(s);
    x = zeros(numel(fn), 1);

    for i = 1:numel(fn)
        % assignin('caller', fn{i}, s.(fn{i}));
        x(i) = s.(fn{i})(ind);
    end

end