% Tim Truster
% 08/25/2013
%
% Template user element routine

% Note: if the user wants to avoid the issues with using global memory
% access (i.e. scripts vs functions), it is recommended to use functions
% that are called within the script, e.g. [lie,nh1,nh3] = e11case1(input).


% % DG Data Load - converts single xl,ul arrays into a left (L) and right 
% % (R) side for separate treatment during calculations. Use only for DG
% % elements.
% if isw ~= 1
% CGtoDGarrays
% 
% inter = elem - (numel - numSI);
% nodeAR = SurfacesI(inter,1);
% nodeBR = SurfacesI(inter,2);
% nodeAL = SurfacesI(inter,3);
% nodeBL = SurfacesI(inter,4);
% 
% nelLP = nelL;
% nelRP = nelR;
% end


%% Global declarations
% Place here any commands that should be executed
% whenever the element routine is called.


%% Task Switch - Perform various element functions
switch isw 
%%    
    case 1 % Setup up elemental degrees of freedom (OPTIONAL)
        
% Purpose: Rearrange order of dofs on element and/or flag certain
%          dofs to be presecribed away
% Called by: pmatin.m
% Notes: Setting values of lie(i,1) zero for i = 1:ndf flags those
%        dofs on all nodes to be unused; global dofs will be
%        prescribed as zero if all elements connected to that node
%        have flagged that dof.
%        Individual dofs on each node are handled by slots lie(i,j)
%        for j=2:nen+1.
%        History variables are allocated using nh1 and nh3.
%        Dof reordering is handled in MatTypeTable in the input
%        file; see NL_Elem8_3d.m for an example.
%        See FEAP pmanual.pdf and pmatin.m for more details.
% Example: NL_Elem8_3d.m
        
        % Example
        if ndf > ndm
            
            for i = ndm+1:ndf
                lie(i,1) = 0;
            end
            
        end

        nh1 = 0; % number of time dependent history variables
        nh3 = 0; % number of time independent history variables
        istv = 0; % number of stresses per node for post-processing
        iste = 0; % number of stresses per element (total number, including all integration points)
        istee = 0; % number of error norm quantities

%%
    case 3 % Stiffness and internal force vector (REQUIRED)
        
        [ElemK,ElemF,hr] = BlankElem03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen);
        

%%        
    case -1 % Boundary tractions (RECOMMENDED)
        
        ElemF = BlankElemm1(mateprop,nodeA,nodeB,elem,edge,traction,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen);


%%
    case 5 % Mass matrix (OPTIONAL)
        
        [ElemM,hr] = BlankElem05(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen);
        
        
%%
    case 6 % Internal force vector (RECOMMENDED)
        
%         [ElemK,ElemF,hr] = BlankElem03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen);
        [ElemF,hr] = BlankElem06(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen);


%%
    case 7 % user boundary tractions (OPTIONAL)

% Purpose: Compute user surface traction for an element; these forces are
%          recomputed at every step (iteration)
% Called by: ploadu.m
% Notes: The ordering of dofs in the external load vector are as follows:
%        ElemF = [node1_dof1 node1_dof2 ... node1_dofndf node2_dof1 ...
%                 node2_dofndf node3_dof1 ... nodenen_dofndf]';
%        Integration of the force can be handled in a standard fashion by
%        rotating the parameteric space of the current element to match a
%        template face over which integration is always performed.
% Example: NL_Elem2_2dM.m        
        

%%
    case 9 %Global error
        
% Purpose: Compute global error residual vector (RHS) for element
% Called by: Implicit_Error_Estimation.m
% Notes: Please see NL_Elem5_2dNSCST.m for an example, where the RHS is
%        computed after the local-implicit error has been calculated.
%        The internal force vector should be the negative of its proper
%        value, e.g. -B'*sigma, due to the definition that the
%        Residual = F_external - F_internal.
%        While nst = ndf*nen in general, typically the assembly process
%        only requires the values for i=1:ndf*nel, where nel is the number
%        of nodes on the current element.
% Example: NL_Elem5_2dNSCST.m
        
        ElemF = zeros(nst,1);
        
        
%%       
    case 11 % Error estimation (OPTIONAL)
        
        ElemE = BlankElem11(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen,numEn);
        
        
%%
    case 12 % Energy (OPTIONAL)
        
        ElemE = BlankElem12(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen);
        
        
%%
    case 15 % Body force calculation (OPTIONAL)
        
        ElemF = BlankElem15(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen);


%%
    case -15 % User Body Force (OPTIONAL)
        
% Purpose: Compute body force for element; these forces are recomputed at 
%          every step (iteration)
% Called by: pbodyfu.m
% Notes: The ordering of dofs in the external load vector are as follows:
%        ElemF = [node1_dof1 node1_dof2 ... node1_dofndf node2_dof1 ...
%                 node2_dofndf node3_dof1 ... nodenen_dofndf]';
%        Body forces are treated as external loads, which can be scaled by
%        a proportional loading parameter.
%        ElemF is initialized as zero inside pbodyf.m and does not need to
%        be reinitialized here.
% Example: NL_Elem2_2dM.m


%%
    case 21 % Stiffness matrix (RECOMMENDED)
        
%         [ElemK,ElemF,hr] = BlankElem03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen);
        [ElemK,hr] = BlankElem21(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen);


%%
    case 24 % Plasticity data (OPTIONAL)
        
% Purpose: Compute plasticity data (for CEE577, single element)
% Called by: FormFE.m
% Notes: Only really tested for a single element mesh.
% Example: Bbar3d_Elem2.m
        
        ElemP = zeros(12,nel);


%%
    case 25 % Stress projection to nodes (RECOMMENDED)
        
        ElemS = BlankElem25(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nelS,nen,npstr);
        

%%        
    case 26 % Element Stress (OPTIONAL)
        
        ElemS = BlankElem26(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen,nestr);
        
        
%%
    case 40 % Initialize history variables (OPTIONAL)
        
        hr = BlankElem40(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen);


%%
    case 51 % Volume stress/strain homogenization (OPTIONAL)
        
% Purpose: Compute volume averaged stresses and strains
% Called by: FormI.m
% Example: NL_Elem8_3d.m (not yet verified)

        
%%
    case 60 % output interface quantities for plotting (OPTIONAL)
        
% Purpose: Compute interface quantities and output; DG elements only
% Called by: FormI.m
% Example: NL_Elem21_2d_2.m


%%
    case 61 % form data structure for interface segments (OPTIONAL)
        
% Purpose: Form data structure for interface segments so that data can be
%          plotted; DG elements only
% Called by: FormIData.m
% Example: NL_Elem21_2d_2.m
        

end %Task Switch
