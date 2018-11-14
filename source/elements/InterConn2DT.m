% Tim Truster
% 08/30/2012
%
% Routine to determine the bounds of integration and other information for
% 2D DG interface elements.
% Generates triangular integration elements for fine scale bubbles.

        xlintL = zeros(2,3);
        xlintR = zeros(2,3);
        
        % Determine bounds of integration, right
        
        if nelR == 4 || nelR == 9
            
            t1 = [(xlR(1,2)-xlR(1,1)); (xlR(2,2)-xlR(2,1)); 0];
            t2 = [(xlR(1,4)-xlR(1,1)); (xlR(2,4)-xlR(2,1)); 0];
            t3 = VecCrossProd(t1,t2);
            aR = VecNormalize(t3);
            hR = sqrt((xlR(1,2)-xlR(1,1))^2+(xlR(2,2)-xlR(2,1))^2)/aR;
            
            drR = 2;
            roR = -1;

            % Upper Limit
            if nodeAR == ElemFlagR(2)
                eR2 = 1;
                xlintR(:,2) = xlR(:,2);
                xlintR(:,3) = xlR(:,3);
            elseif nodeAR == -1 %no enrichment but DG instead
                xy = xlL(:,1);
                POUxi = POU_Coord(xy(1),xy(2),xlR,1,nelR);
                eR2 = POUxi(1);
                xlintR(:,2) = xy;
            elseif nelR == 9 && nodeAR == ElemFlagR(5)
                eR2 = 0;
            end
            % Lower Limit
            if nodeBR == ElemFlagR(1)
                eR1 = -1;
                xlintR(:,1) = xlR(:,1);
                xlintR(:,3) = xlR(:,4);
            elseif nodeBR == -1 %no enrichment but DG instead
                xy = xlL(:,2);
                POUxi = POU_Coord(xy(1),xy(2),xlR,1,nelR);
                eR1 = POUxi(1);
                xlintR(:,1) = xy;
                xlintR(:,3) = xlR(:,3);
            elseif nelR == 9 && nodeBR == ElemFlagR(5)
                eR1 = 0;
            end
        
        elseif nelR == 3 || nelR == 6
            
            t1 = [(xlR(1,2)-xlR(1,1)); (xlR(2,2)-xlR(2,1)); 0];
            t2 = [(xlR(1,3)-xlR(1,1)); (xlR(2,3)-xlR(2,1)); 0];
            t3 = VecCrossProd(t1,t2);
            aR = VecNormalize(t3);
            aR = aR/2;
            hR = sqrt((xlR(1,2)-xlR(1,1))^2+(xlR(2,2)-xlR(2,1))^2)/aR;
            
            drR = 1;
            roR = 0;

            % Upper Limit
            if nodeAR == ElemFlagR(2)
                eR2 = 1;
                xlintR(:,2) = xlR(:,2);
            elseif nodeAR == -1 %no enrichment but DG instead
                xy = xlL(:,1);
                POUxi = POU_Coord(xy(1),xy(2),xlR,1,nelR);
                eR2 = POUxi(1);
                xlintR(:,2) = xy;
            elseif nelR == 6 && nodeAR == ElemFlagR(4)
                eR2 = 1/2;
            end
            % Lower Limit
            if nodeBR == ElemFlagR(1)
                eR1 = 0;
                xlintR(:,1) = xlR(:,1);
            elseif nodeBR == -1 %no enrichment but DG instead
                xy = xlL(:,2);
                POUxi = POU_Coord(xy(1),xy(2),xlR,1,nelR);
                eR1 = POUxi(1);
                xlintR(:,1) = xy;
            elseif nelR == 6 && nodeBR == ElemFlagR(4)
                eR1 = 1/2;
            end
            xlintR(:,3) = xlR(:,3);
        
        end
        
        % Determine bounds of integration, left
        
        if nelL == 4 || nelL == 9
            
            t1 = [(xlL(1,2)-xlL(1,1)); (xlL(2,2)-xlL(2,1)); 0];
            t2 = [(xlL(1,4)-xlL(1,1)); (xlL(2,4)-xlL(2,1)); 0];
            t3 = VecCrossProd(t1,t2);
            aL = VecNormalize(t3);
            hL = sqrt((xlL(1,2)-xlL(1,1))^2+(xlL(2,2)-xlL(2,1))^2)/aL;
            
            drL = 2;
            roL = -1;

            % Upper Limit
            if nodeAL == ElemFlagL(1)
                eL1 = -1;
                xlintL(:,1) = xlL(:,1);
                xlintL(:,3) = xlL(:,4);
            elseif nodeAL == -1 %no enrichment but DG instead
                xy = xlR(:,2);
                POUxi = POU_Coord(xy(1),xy(2),xlL,1,nelL);
                eL1 = POUxi(1);
                xlintL(:,1) = xlL(:,1);
            elseif nelL == 9 && nodeAL == ElemFlagL(5)
                eL1 = 0;
            end
            % Lower Limit
            if nodeBL == ElemFlagL(2)
                eL2 = 1;
                xlintL(:,2) = xlL(:,2);
                xlintL(:,3) = xlL(:,3);
            elseif nodeBL == -1 %no enrichment but DG instead
                xy = xlR(:,1);
                POUxi = POU_Coord(xy(1),xy(2),xlL,1,nelL);
                eL2 = POUxi(1);
                xlintL(:,2) = xlL(:,2);
                xlintL(:,3) = xlL(:,4);
            elseif nelL == 9 && nodeB == ElemFlagL(5)
                eL2 = 0;
            end
        
        elseif nelL == 3 || nelL == 6
            
            t1 = [(xlL(1,2)-xlL(1,1)); (xlL(2,2)-xlL(2,1)); 0];
            t2 = [(xlL(1,3)-xlL(1,1)); (xlL(2,3)-xlL(2,1)); 0];
            t3 = VecCrossProd(t1,t2);
            aL = VecNormalize(t3);
            aL = aL/2;
            hL = sqrt((xlL(1,2)-xlL(1,1))^2+(xlL(2,2)-xlL(2,1))^2)/aL;
            
            drL = 1;
            roL = 0;

            % Upper Limit
            if nodeAL == ElemFlagL(1)
                eL1 = 0;
                xlintL(:,1) = xlL(:,1);
            elseif nodeAL == -1 %no enrichment but DG instead
                xy = xlR(:,2);
                POUxi = POU_Coord(xy(1),xy(2),xlL,1,nelL);
                eL1 = POUxi(1);
                xlintL(:,1) = xy;  
            elseif nelL == 6 && nodeAL == ElemFlagL(4)
                eL1 = 1/2;
            end
            % Lower Limit
            if nodeBL == ElemFlagL(2)
                eL2 = 1;
                xlintL(:,2) = xlL(:,2);
            elseif nodeBL == -1 %no enrichment but DG instead
                xy = xlR(:,1);
                POUxi = POU_Coord(xy(1),xy(2),xlL,1,nelL);
                eL2 = POUxi(1);
                xlintL(:,2) = xy;  
            elseif nelL == 6 && nodeBL == ElemFlagL(4)
                eL2 = 1/2;
            end
            xlintL(:,3) = xlL(:,3);  
        
        end