function dE = dHpcs(g, c)
% gets the attraction energy differential for cell c.
% E(c) = d where d is the distance from cell c to the line on PCs.
dE = zeros(2*length(g.verts),1);

if g.is_LImodel

absolute_top_position = 0; % indicates if the position to pull to is absolute or depends on the PCs current position
if ismember(c, g.LImodel.high_delta_cells)
    
   vidx = getVerts(g,c);
   vn = length(vidx); % number of vertices in cell c
   vert = getRelativePosition(g,vidx);
   cpos = g.centroid(c,:);
   A = cellarea(g, c);
   xc = cpos(1); yc = cpos(2);
   dist_thresh = 2*sqrt(g.A0/pi);
   
   if absolute_top_position
       y_top = 0;
   else
       y_top = mean(g.centroid(g.top_boundary_cells, 2));
   end
   D = yc - y_top; % the distance from top boundary
   if abs(D) < dist_thresh, return; end

   for k=1:vn
       kn = mod(k,vn)+1; % next vertex
       kp = mod(k-2,vn)+1; % previous vertex
       r0 = vert(kp,:); r1 = vert(k,:); r2 = vert(kn,:); % positions of vertices k-1, k, k+1 respectively
       xp = r0(1); x = r1(1); xn = r2(1); yp = r0(2); y = r1(2); yn = r2(2);
       dAdx = 0.5*( yn-yp );
       dAdy = 0.5*( xp-xn );
       dYcmdx = ( -dAdx*yc + (1/6)*(yn*(y+yn)-yp*(y+yp)) )/A;
       dYcmdy = ( -dAdy*yc + (1/6)*(2*y*(xp-xn)+x*(yn-yp)+yp*xp-yn*xn) )/A;
       dE(2*vidx(k)-1) = dE(2*vidx(k)-1) - 2*D*dYcmdx;
       dE(2*vidx(k)) = dE(2*vidx(k)) - 2*D*dYcmdy;
   end
 
end
dE = g.paras(5)*dE/(g.A0);

end
end