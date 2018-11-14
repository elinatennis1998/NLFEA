% 06/16/2013
% Tim Truster
% Implicit Error Calculation, segregated from NL_FEA_Program

        if implicon == 1
            
            if subm == 1

        %     Variables to Pass for array dimension declarations

                minc = 2; % Number of Sub-Cells in each direction

        %         celn = (minc+1)^2; %(minc+1)*(minc+1)! Number of Nodes in Patch
                celn = (2*minc+1)^2; %(minc+1)*(minc+1)! Number of Nodes in Patch

                nstar = 2; %(m+1)*(m+1)*ndf	% Number of DOFs in Patch

                cel = minc^2; %(m-1)**2			% Interior nodes in Patch

                maxel = 8;
                ndfs = 2;
                strwek = residflag;
                strong = strwek;
                x = Coordinates';
                NodesOnElement = NodesOnElement';
                nen1 = nen+1;
                u = squeeze(DispList(:,:,step));
                d = MateT';
                IQTY = 1;
                rho = 1;
                [i,nummat] = size(d);
                d = [d; rho*ones(1,nummat); IQTY*ones(1,nummat); iprob*ones(1,nummat)];
                if minc > 0
                    submesh
                end
                NodesOnElement = NodesOnElement';

            end

        else

            minc = 2; % Number of Sub-Cells in each direction

        %         celn = (minc+1)^2; %(minc+1)*(minc+1)! Number of Nodes in Patch
            celn = (2*minc+1)^2; %(minc+1)*(minc+1)! Number of Nodes in Patch

            nstar = 2; %(m+1)*(m+1)*ndf	% Number of DOFs in Patch

            cel = minc^2; %(m-1)**2			% Interior nodes in Patch

            maxel = 8;
            ndfs = 2;
            strwek = residflag;
            strong = strwek;
            x = Coordinates';
            NodesOnElement = NodesOnElement';
            nen1 = nen+1;
            gensubmeshorig
            NodesOnElement = NodesOnElement';

        end

%     if subm == 1
        
%         Pd = zeros(neq, 1);
%         isw = 7;
%         FormFE
%         
        %Assemble Stiffness Routine

        Fd = zeros(neq, 1);
        isw = 9;
        FormFE
%         isw = 21;
%         FormFE

        %Solve Matrix System for FE Solution

        qbar2 = Kdd11\Fd1;

        u = zeros(numnp,ndf);

        for node = 1:numnp
            for dir = 1:ndf
                gDOF = NDOFT(node, dir);
                if gDOF <= neq
                    u(node, dir) = qbar2(gDOF,1);
                else
%                     u(node, dir) = gBC(gDOF - neq);
                end
            end
        end

        u = u';

        projectCtoF

%         plotModelContP(xs', up(1,:)', ixs', numels, nen, 12, 1, 1, 'Global Error - e_T')
%         plotModelContP(Coordinates, u(1,:)', ix, numel, nen, 13, 1, 1, 'Global Error - e_T')
%         plotModelErroP(Coordinates, Node_U_V(:,1), ix, numel, nen, 2, 1, 1, 'Standard Error - u_x',iprob,1)
        
        implicitnormsT
        if implicon == 1
        Ieffvals2 = sqrt(ufieff./Ieffvals(:,3));
        end

        uo = up;        
        u = squeeze(DispList(:,:,step));

        projectCtoF
        
% Stored variables at this point:
% u = Node_U_V
% ut = total error estimate: sum of fine scale and coarse scale errors
% up = Node_U_V projected onto submesh
% uo = global errors projected onto submesh
% xs = refined submesh nodal coordinates

% % Plot of errors on updated deformed mesh
%  plotModelCont(xs(:,1:numnps)'+ut(1:2,:)'+up(1:2,:)', ut(1,:)', ixs', numels, nen, 1, 1, 1, '')
%  plotModelCont(xs(:,1:numnps)'+ut(1:2,:)'+up(1:2,:)', ut(2,:)', ixs', numels, nen, 2, 1, 1, '')
%  plotModelCont(xs(:,1:numnps)'+ut(1:2,:)'+up(1:2,:)', ut(3,:)', ixs', numels, nen, 3, 1, 1, '')
% % Plot of bubble fine-scale H1 error norms
% plotIeffContP(Coordinates, Ieffvals(:,2), ix, numel, nen, 4, 1, 1, 'Local Effectivity Indices')
% % Plot of submesh fine-scale H1 error norms
% plotIeffContP(Coordinates, ufieff', ix, numel, nen, 5, 1, 1, 'Local Effectivity Indices')
% % Plot of total H1 error norms
% plotIeffContP(Coordinates, utieff', ix, numel, nen, 8, 1, 1, 'Local Effectivity Indices')
% % Plot of bubble fine scale errors on deformed mesh
% plotModelContBP(Coordinates+Node_U_V(:,1:2),0*u', bubblevals(:,1), ix, numel, nen, 6, 1, 1, 'Global Error - e_T')
% plotModelContBP(Coordinates+Node_U_V(:,1:2),0*u', bubblevals(:,2), ix, numel, nen, 7, 1, 1, 'Global Error - e_T')