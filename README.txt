Upon opening MATLAB, the user should load the 'source' directory into the path:
addpath(genpath(['full_path_to_repository' NLFEA\source']))

Revised input variable names 9/10/17 to match DEIP:
NLFEA_old_name              DEIP_new_name
NodeTable                   Coordinates
ix(:,1:nen)                 NodesOnElement
ix(:,nen1)                  RegionOnElement