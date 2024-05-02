function [in_central_name, out_central_name] = In_out_degree_centrality(W,name)
    % Sum of the columns is the in-degree
    in_w = sum(W);
    % Sum of the rows is the out-degree
    out_w = sum(W');

    % Find the three most central sectors
    max_3_in_w = Find_idx(in_w,3);
    max_3_out_w = Find_idx(out_w,3);

    % Find the name
    in_central_name = [name(max_3_in_w(1)), name(max_3_in_w(2)), name(max_3_in_w(3))];
    out_central_name = [name(max_3_out_w(1)), name(max_3_out_w(2)), name(max_3_out_w(3))];
end