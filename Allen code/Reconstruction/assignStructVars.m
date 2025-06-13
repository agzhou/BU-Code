% assign fields from a structure into the workspace that calls this function
% adapted from https://stackoverflow.com/a/29958475

function assignStructVars(s)

    fn = fieldnames(s);
    for i = 1:numel(fn)
        assignin('caller', fn{i}, s.(fn{i}));
    end

end