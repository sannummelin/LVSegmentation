function normals = compute_polygon_normals(x, y)
    % input: x, y - column vectors of polygon vertex coordinates (N x 1)
    % output:normals - (N x 2) matrix of normal vectors at each vertex

    % ensure column vectors
    x = x(:);
    y = y(:);

    % close the polygon if not already
    if ~(x(1) == x(end) && y(1) == y(end))
        x(end+1) = x(1);
        y(end+1) = y(1);
    end

    % compute edge vectors
    dx = diff(x);
    dy = diff(y);

    % compute normal vectors to the edges by rotating 90 degrees
    nx = dy;
    ny = -dx;

    % normalize edge normals to unit length
    L = sqrt(nx.^2 + ny.^2);
    nx = nx ./ L;
    ny = ny ./ L;

    % duplicate the first normal at the end for averaging at vertices
    nx(end+1) = nx(1);
    ny(end+1) = ny(1);

    % number of vertices in the closed polygon
    N = length(x);

    % initialize output normal vector array
    normals = zeros(N, 2);

    % interpolate normals to each vertex by averaging adjacent edge normals
    for i = 1:N
        % get previous and current edge normals
        n1 = [nx(mod(i-2, N)+1), ny(mod(i-2, N)+1)];
        n2 = [nx(i), ny(i)];

        % average the two normals
        n_avg = (n1 + n2) / 2;

        % normalize the averaged normal vector
        n_avg = n_avg / norm(n_avg);

        % store the result
        normals(i, :) = n_avg;
    end
end