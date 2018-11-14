function [ tetra_node2, tetra_num2, node_num2 ] = tet_mesh_q2l ( tetra_num1, tetra_node1, node_num1 )

%*****************************************************************************80
%
%% TET_MESH_Q2L is the main program for TET_MESH_Q2L.
%
%  Discussion:
%
%    TET_MESH_Q2L converts a quadratic tet mesh to a linear one.
%
%    A quadratic tet mesh is assumed to consist of 10-node
%    tetrahedrons.  This routine rearranges information so as to 
%    define a 4-node tet mesh.
%
%  Usage:
%
%    tet_mesh_q2l ( 'prefix' )
%
%    where
%
%    * 'prefix'_nodes.txt contains nodal coordinates (not used by this function);
%    * 'prefix'_elements.txt contains the element definitions;
%    * 'prefix'_q2l_elements.txt contains the new element definitions;
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    17 August 2009
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string PREFIX, the common filename prefix.
%

%
%  Set the number of linear tetrahedrons:
%

  [ node_num2, tetra_num2 ] = tet_mesh_order10_to_order4_size ( ...
    node_num1, tetra_num1 );

%
%  Convert the data.
%
  tetra_node2 = tet_mesh_order10_to_order4_compute ( tetra_num1, ...
    tetra_node1, tetra_num2 );

end
function tetra_node2 = tet_mesh_order10_to_order4_compute ( tetra_num1, ...
  tetra_node1, tetra_num2 )

%*****************************************************************************80
%
%% TET_MESH_ORDER10_TO_ORDER4_COMPUTE linearizes a quadratic tet mesh.
%
%  Discussion:
%
%    A quadratic tet mesh is assumed to consist of 10-node
%    tetrahedrons.
%
%    This routine rearranges the information so as to define a 4-node
%    tet mesh.
%
%    The same nodes are used, but there are 8 times as many
%    tetrahedrons.
%
%    The node ordering for the quadratic tetrahedron is somewhat
%    arbitrary.  In the current scheme, the vertices are listed
%    first, followed by the 6 midside nodes.  Each midside node
%    may be identified by the two vertices that bracket it.  Thus,
%    the node ordering may be suggested by:
%
%      1  2  3  4 (1+2) (1+3) (1+4) (2+3) (2+4) (3+4)
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    22 January 2007
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Anwei Liu, Barry Joe,
%    Quality Local Refinement of Tetrahedral Meshes Based
%    on 8-Subtetrahedron Subdivision,
%    Mathematics of Computation,
%    Volume 65, Number 215, July 1996, pages 1183-1200.
%
%  Parameters:
%
%    Input, integer TETRA_NUM1, the number of tetrahedrons in the quadratic
%    tet mesh.
%
%    Input, integer TETRA_NODE1(10,TETRA_NUM1), the quadratic tet mesh.
%
%    Input, integer TETRA_NUM2, the number of tetrahedrons in the linear
%    tet mesh.  TETRA_NUM2 = 8 * TETRA_NUM1.
%
%    Output, integer TETRA_NODE2(4,TETRA_NUM2), the linear tet mesh.
%
  tetra_node2 = zeros ( 4, tetra_num2 );
  
  tetra2 = 0;

  for tetra1 = 1 : tetra_num1

    n1 = tetra_node1(1,tetra1);
    n2 = tetra_node1(2,tetra1);
    n3 = tetra_node1(3,tetra1);
    n4 = tetra_node1(4,tetra1);
    n5 = tetra_node1(5,tetra1);
    n6 = tetra_node1(6,tetra1);
    n7 = tetra_node1(7,tetra1);
    n8 = tetra_node1(8,tetra1);
    n9 = tetra_node1(9,tetra1);
    nx = tetra_node1(10,tetra1);

%     tetra2 = tetra2 + 1;
%     tetra_node2(1:4,tetra2 ) = [ n1, n5, n6, n7 ]';
%     tetra2 = tetra2 + 1;
%     tetra_node2(1:4,tetra2 ) = [ n2, n5, n8, n9 ]';
%     tetra2 = tetra2 + 1;
%     tetra_node2(1:4,tetra2 ) = [ n3, n6, n8, n9 ]';
%     tetra2 = tetra2 + 1;
%     tetra_node2(1:4,tetra2 ) = [ n4, n7, n9, nx ]';
%     tetra2 = tetra2 + 1;
%     tetra_node2(1:4,tetra2)  = [ n5, n6, n7, n9 ]';
%     tetra2 = tetra2 + 1;
%     tetra_node2(1:4,tetra2 ) = [ n5, n6, n8, n9 ]';
%     tetra2 = tetra2 + 1;
%     tetra_node2(1:4,tetra2 ) = [ n6, n7, n9, nx ]';
%     tetra2 = tetra2 + 1;
%     tetra_node2(1:4,tetra2 ) = [ n6, n8, n9, nx ]';
% % TJT
%     tetra2 = tetra2 + 1;
%     tetra_node2(1:4,tetra2 ) = [ n1, n5, n7, n8 ]';
%     tetra2 = tetra2 + 1;
%     tetra_node2(1:4,tetra2 ) = [ n2, n5, n9, n6 ]';
%     tetra2 = tetra2 + 1;
%     tetra_node2(1:4,tetra2 ) = [ n3, n6, n7, nx ]';
%     tetra2 = tetra2 + 1;
%     tetra_node2(1:4,tetra2 ) = [ n4, n8, nx, n9 ]';
%     tetra2 = tetra2 + 1;
%     tetra_node2(1:4,tetra2)  = [ n5, n6, n7, n8 ]';
%     tetra2 = tetra2 + 1;
%     tetra_node2(1:4,tetra2 ) = [ n5, n6, n8, n9 ]';
%     tetra2 = tetra2 + 1;
%     tetra_node2(1:4,tetra2 ) = [ n6, n7, n8, nx ]';
%     tetra2 = tetra2 + 1;
%     tetra_node2(1:4,tetra2 ) = [ n6, n8, n9, nx ]';
    tetra2 = tetra2 + 1;
    tetra_node2(1:4,tetra2 ) = [ n1, n5, n7, n8 ]';
    tetra2 = tetra2 + 1;
    tetra_node2(1:4,tetra2 ) = [ n2, n5, n9, n6 ]';
    tetra2 = tetra2 + 1;
    tetra_node2(1:4,tetra2 ) = [ n3, n7, n6, nx ]';
    tetra2 = tetra2 + 1;
    tetra_node2(1:4,tetra2 ) = [ n4, n8, nx, n9 ]';
    tetra2 = tetra2 + 1;
    tetra_node2(1:4,tetra2)  = [ n5, n7, n8, n9 ]';
    tetra2 = tetra2 + 1;
    tetra_node2(1:4,tetra2 ) = [ n5, n7, n9, n6 ]';
    tetra2 = tetra2 + 1;
    tetra_node2(1:4,tetra2 ) = [ n7, n8, n9, nx ]';
    tetra2 = tetra2 + 1;
    tetra_node2(1:4,tetra2 ) = [ n7, n6, nx, n9 ]';

  end

  return
end
function [ node_num2, tetra_num2 ] = tet_mesh_order10_to_order4_size ( ...
  node_num1, tetra_num1 )

%*****************************************************************************80
%
%% TET_MESH_ORDER10_TO_ORDER4_SIZE sizes a linear tet mesh from a quadratic one.
%
%  Discussion:
%
%    A linear (4 node) tet mesh can be derived from a quadratic
%    (10 node) tet mesh using the same set of nodes, but reassigning
%    the nodes of each quadratic tet among 8 linear subtets.
%
%    This routine returns the number of nodes and tetrahedra in the
%    linear mesh.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    22 January 2007
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Anwei Liu, Barry Joe,
%    Quality Local Refinement of Tetrahedral Meshes Based
%    on 8-Subtetrahedron Subdivision,
%    Mathematics of Computation,
%    Volume 65, Number 215, July 1996, pages 1183-1200.
%
%  Parameters:
%
%    Input, integer TETRA_NUM1, the number of tetrahedrons in the
%    quadratic mesh.
%
%    Input, integer NODE_NUM1, the number of nodes in the quadratic mesh.
%
%    Output, integer TETRA_NUM2, the number of tetrahedrons in the
%    linear mesh.
%
%    Output, integer NODE_NUM2, the number of nodes for the linear mesh.
%
  node_num2 = node_num1;
  tetra_num2 = 8 * tetra_num1;

  return
end
