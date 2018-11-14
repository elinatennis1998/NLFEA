% 4/29/2015
%
% Sample for: L_Elem1_2d
% Notes:
%
% Input File: Two Trianglular Elements Under Axial Load
%
% This input file should be run prior to executing the FEA_Program routine.
%
% Format of required input:
%
%   numnp:           = number of nodes in the mesh (length(Coordinates))
%
%   numel:           = number of elements in the mesh
%
%   nummat:
%
%   ndm:
%
%   ndf:
%
%   nen:             = maximum number of nodes per element (4)
%
%   nen1:
%
%   PSPS:            = flag for plane stress ('s') or plane strain ('n')
%
%   Coordinates:       = table of mesh nodal coordinates defining the
%                      geometry of the mesh; format of the table is as
%                      follows:
%                          Nodes  |             x-coord  y-coord
%                          n1     |  Coordinates = [x1     y1
%                          n2     |               x2     y2
%                          ...    |               ..     ..
%                          nnumnp |               xnumnp ynumnp];
%
%   ix:              = table of mesh connectivity information, specifying
%                      how nodes are attached to elements and how materials
%                      are assigned to elements; entries in the first nen
%                      columns correspond to the rows of Coordinates
%                      representing the nodes attached to element e;
%                      entries in the last nen+1 column are rows from MateT
%                      signifying the material properties assigned to
%                      element e; format of the table is as follows:
%                          Elements  |         n1    n2    n3    n4   mat
%                          e1        |  ix = [e1n1  e1n2  e1n3  e1n4 e1mat
%                          e2        |        e2n1  e2n2  e2n3  e2n4 e2mat
%                          ...       |         ..    ..    ..    ..   ..
%                          enumel    |        values for element numel   ];
%
%   MateT:           = table of mesh material properties for each distinct
%                      set of material properties; these sets are
%                      referenced by element e by setting the value of
%                      ix(e,nen+1) to the row number of the desired
%                      material set; format of the table is as follows:
%                          Materials  |           E   v   t
%                          mat1       |  MateT = [E1  v1  t1
%                          mat2       |           E2  v2  t2
%                          ...        |           ..  ..  ..];
%
%   MatTypeTable: = [a b c ]
%                   b element subroutine number(21 for nonlinear DG interface)
%                   c is the flag for 1 nonlinear or 0 linear
%   AlgoType:
%
%   OptFlag:
%
%   BCLIndex:        = list of the number of boundary conditions and loads
%                      applied to the mesh; first entry is the number of
%                      prescribed displacements at nodes; second entry is
%                      the number of nodal forces
%
%   NodeBC:          = table of prescribed nodal displacement boundary
%                      conditions; it contains lists of nodes, the
%                      direction of the displacement prescribed (x=1, y=2),
%                      and the value of the displacement (set 0 for fixed
%                      boundary); the length of the table must match the
%                      entry in BCLIndex(1), otherwise an error will result
%                      if too few conditions are given or extra BCs will be
%                      ignored in the model input module;  format of the 
%                      table is as follows:
%                          BCs  |            nodes direction value
%                          bc1  |  NodeBC = [bc1n   bc1dir   bc1u
%                          bc2  |            bc2n   bc2dir   bc2u
%                          ...  |             ..     ..       .. ];
%
%   NodeLoad:        = table of prescribed nodal forces; it contains lists 
%                      of nodes, the direction of the force prescribed 
%                      (x=1, y=2), and the value of the force; the length 
%                      of the table must match the entry in BCLIndex(2), 
%                      otherwise an error will result if too few conditions
%                      are given or extra loads will be ignored in the 
%                      model input module; format of the table is as
%                      follows:
%                          Loads  |              nodes direction value
%                          P1     |  NodeLoad = [ P1n    P1dir    P1P
%                          P2     |               P2n    P2dir    P2P
%                          ...    |               ..     ..       .. ];
%
%   iprob: = 6 for body force problem
%
%   numSI:
%
% The following numbering convention is used for 4-node quadrilateral
% elements:
%
%           4 -------------- 3
%           |                |
%           |                |
%           |                |
%           |                |
%           |                |
%           |                |
%           1 -------------- 2
%

% Arbitrary data for assistance in defining the mesh
L = 40;
H = 20;
q = 0.2;
P = L/2*q;

% Mesh Nodal Coordinates
Coordinates = [0 0
             L 0
             0 H
             L H];
numnp = length(Coordinates);

% Mesh Element Connectivities
NodesOnElement = [1 4 3 0
      1 2 4 0];
RegionOnElement = [1 1]';
nen = 4;
nen1 = nen + 1;
numel = 2;
iprob = 0;

% Mesh Boundary Conditions and Loads
NodeBC = [1 2 0
          2 1 0
          2 2 0];
NodeLoad = [3 2 1000*P
            4 2 1000*P];
numBC = 3;
numNodalF = 2;
SurfacesL = [ 3 4 1 2    0     -2*P/L     0];
% numSL = 1;
numSL = 0;

% Mesh Material Properties
AlgoType = [-1; 1; 0];
OptFlag = [0 1 1 0 0 0 0]';
OptFlag = [0 1 1 1 1 1 1]';
young = 1e6;
pois = .25;
thick = 1;
MateT = [young pois thick];
MatTypeTable = [1; 1; 0];

nummat = 1;
numSI = 0;
ProbType = [numnp numel nummat 2 2 nen];
