dim = 2;
level = 6;
order = 0; % not used by Clenshaw-Curtis rule
[weights, points] = tsgMakeQuadrature(dim, 'clenshaw-curtis', 'level', level, order);

I = weights' * (exp(-points(:,1).^2) .* cos(points(:,2))); % I is the approximated quadrature
E = 2.513723354063905e+00; % E is the "exact" solution computed to 16 decimal places
E = abs(I - E);
%E = abs(I - quad(@(x)(exp(-x.^2)), -1, 1, 1.E-14)*quad(@(x)(cos(x)), -1, 1, 1.E-14));

disp(['----------------------------------------------------------------------------']);
disp([' Example 1:  integrate f(x,y) = exp(-x^2) * cos(y), using clenshaw-curtis level nodes']);
disp(['    at level ',num2str(level)]);
disp(['    grid has ',num2str(size(points,1)),' points']);
disp(['    integral: ',num2str(I,16)]);
disp(['    error: ',num2str(E,16)]);
disp([' ']);

level = 7;
[weights, points] = tsgMakeQuadrature(dim, 'clenshaw-curtis', 'level', level, order);

I = weights' * (exp(-points(:,1).^2) .* cos(points(:,2))); % I is the approximated quadrature
E = 2.513723354063905e+00; % E is the "exact" solution computed to 16 decimal places
E = abs(I - E);

disp(['    at level ',num2str(level)]);
disp(['    grid has ',num2str(size(points,1)),' points']);
disp(['    integral: ',num2str(I,16)]);
disp(['    error: ',num2str(E,16),' (rounded to 14 decimal places)']);
disp([' ']);