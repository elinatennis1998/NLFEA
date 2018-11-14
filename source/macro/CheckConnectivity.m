function [flag,inc21,inc22,start1,start2] = CheckConnectivity(Corners1,Corners2,minc,inc21,inc22)
% [flag,inc21,inc22,start1,start2,corn1,corn2,trans] =
% CheckConnectivity(FacePoints1,U1,V1,FacePoints2,U2,V2,n,p,inc21,inc22)
%
% Tim Truster
% CEE Graduate Student
% UIUC
% 1/16/2009
%
%  This subroutine compares two NURBS patch faces to ensure that they are
%  compatible. Compatibility is determined by ensuring that all node
%  weighted-control point coordinates are the same along the face and that
%  the knot vectors parameterizing each face are simply linear combinations
%  of each other. The two faces can be separated such that a linear
%  translation would make them compatible. The proper orientation of the
%  respective parameterizations and this linear translation are determined
%  by the subroutine.
%
%  Note: flag is returned as 2 if face2 needs translated to match face1

% %Extract corner nodes from face
% Corners1 = zeros(4,3);
% Corners2 = zeros(4,3);
% n = (minc +1)*ones(2,3);
% 
% for i = 1:3
% 	Corners1(1,i) = FacePoints1(1,1,i);
% 	Corners1(2,i) = FacePoints1(1,n(1,2),i);
% 	Corners1(3,i) = FacePoints1(n(1,1),n(1,2),i);
% 	Corners1(4,i) = FacePoints1(n(1,1),1,i);
% 	Corners2(1,i) = FacePoints2(1,1,i);
% 	Corners2(2,i) = FacePoints2(1,n(2,2),i);
% 	Corners2(3,i) = FacePoints2(n(2,1),n(2,2),i);
% 	Corners2(4,i) = FacePoints2(n(2,1),1,i);
% end

%Check whether any corners on Face2 match corner1 on Face1
flag = 0;
i = 0;

while i < 4 && flag == 0
	i = i + 1;
    for j = 1:3
		if abs(Corners1(1,j) - Corners2(i,j)) < 10e-9
			flag = flag + 1;
		end
    end
end

if flag == 3 %one node matches
	Match1 = i;
	flag = 0;
    %move to next counter-clockwise (ccw) corner
    if Match1 == 1
        Match12 = 4;
    else
        Match12 = Match1 - 1;
    end
	for j = 1:3
		if abs(Corners1(2,j) - Corners2(Match12,j)) < 10e-9
			flag = flag + 1;
		end
	end
	if flag == 3 %second ccw corner matched
			Match2 = Match12;
			Match22 = -1;
    else %move to next clockwise (cw) corner
        flag = 0;
		if Match1 == 4
			Match12 = 1;
		else
			Match12 = Match1 + 1;
		end
		for j = 1:3
			if abs(Corners1(2,j) - Corners2(Match12,j)) < 10e-9
				flag = flag + 1;
			end
		end
		if flag == 3 %second cw corner matched
				Match2 = Match12;
				Match22 = 1;
		end
	end
	if flag == 3 %second ccw or cw corner matched
		node = Match2 + Match22; %move to next corner
		if node > 4 || node < 1 %check if corner needs looped back around
			node = node -4*Match22;
		end
		flag = 0;
		for j = 1:3
			if abs(Corners1(3,j) - Corners2(node,j)) < 10e-9
				flag = flag + 1;
			end
		end
		if flag == 3 %third corner matched
		end
	end
	if flag == 3 %third ccw or cw corner matched
		node = node + Match22; %move to next corner
		if node > 4 || node < 1 %check if corner needs looped back around
			node = node -4*Match22;
		end
		flag = 0;
		for j = 1:3
			if abs(Corners1(4,j) - Corners2(node,j)) < 10e-9
				flag = flag + 1;
			end
		end
		if flag == 3 %fourth corner matched
		end
		if flag < 3 % 3 coplanar points and one not coplanar
            flag = 0;
            start1 = 0;
            start2 = 0;
			return
		end
	end
end

%Assign corner flags, loop increments, and starting indices based on the
%orientation of the corners

%Match1 identifies the corner on Face2 matching corner1 on Face1
%Match2 identifies the corner on Face2 matching corner2 on Face1

%corn1 identifies which corner on Face1 matches corner1 on Face2
%corn2 identifies which corner on Face2 matches corner1 on Face1

%Note: ALL corner references are created w.r.t. Face(1) serving as the
%master face, which is oriented as follows:
%    4 - 3
%    | 1 |
%    1 - 2
%So, to see the view of Face1 from Face2, simply take the orientation map
%of Face2 to Face1 and rotate/flip it to align it with the master above;
%then perform THE SAME reorientation to Face1, and you will get the map for
%Face1 (i.e. invert the map)

swap = 0;

if Match1 == 1 %origin of parameterization is the same for both faces
	start1 = 1; %     Face2Map       Face2   to   Face1
	start2 = 1; %       4 - 3        4 - 3        4 - 3
	inc1 = 1;   %       | 2 |        | 2 |        | 1 |
	inc2 = 1;   %       1 - 2        1 - 2        1 - 2
	if Match2 == 4 %parameterization is left-handed
%         for j = 1:3
%             FacePoints2(:,:,j) = FacePoints2(:,:,j)';
%         end
		swap = 1;    %Face2Map       Face2   to   Face1
        temp = inc2; %  2 - 3        4 - 3        2 - 3
        inc2 = inc1; %  | 2 |        | 2 |        | 1 |
        inc1 = temp; %  1 - 4        1 - 2        1 - 4
        temp = start2;
        start2 = start1;
        start1 = temp;
	end
elseif Match1 == 2
	start1 = 1;      %Face2Map       Face2   to   Face1
	start2 = minc+1; %  3 - 4        4 - 3        3 - 4
	inc1 = 1;   %       | 2 |        | 2 |        | 1 |
	inc2 = -1;  %       2 - 1        1 - 2        2 - 1
	if Match2 == 3 %parameterization is left-handed
%         for j = 1:3
%             FacePoints2(:,:,j) = FacePoints2(:,:,j)';
%         end
		swap = 1;    %Face2Map       Face2   to   Face1
        temp = inc2; %  1 - 4        4 - 3        3 - 2
        inc2 = inc1; %  | 2 |        | 2 |        | 1 |
        inc1 = temp; %  2 - 3        1 - 2        4 - 1
        temp = start2;
        start2 = start1;
        start1 = temp;
	end
elseif Match1 == 3 %origin of parameterization is katty-corner on the faces
	start1 = minc+1; %Face2Map       Face2   to   Face1
	start2 = minc+1; %  2 - 1        4 - 3        2 - 1
	inc1 = -1;  %       | 2 |        | 2 |        | 1 |
	inc2 = -1;  %       3 - 4        1 - 2        3 - 4
	if Match2 == 2 %parameterization is left-handed
%         for j = 1:3
%             FacePoints2(:,:,j) = FacePoints2(:,:,j)';
%         end
		swap = 1;    %Face2Map       Face2   to   Face1
        temp = inc2; %  4 - 1        4 - 3        4 - 1
        inc2 = inc1; %  | 2 |        | 2 |        | 1 |
        inc1 = temp; %  3 - 2        1 - 2        3 - 2
        temp = start2;
        start2 = start1;
        start1 = temp;
	end
else
	start1 = minc+1; %Face2Map       Face2   to   Face1
	start2 = 1; %       1 - 2        4 - 3        1 - 2
	inc1 = -1;  %       | 2 |        | 2 |        | 1 |
	inc2 = 1;   %       4 - 3        1 - 2        4 - 3
	if Match2 == 1 %parameterization is left-handed
%         for j = 1:3
%             FacePoints2(:,:,j) = FacePoints2(:,:,j)';
%         end
		swap = 1;    %Face2Map       Face2   to   Face1
        temp = inc2; %  3 - 2        4 - 3        1 - 4
        inc2 = inc1; %  | 2 |        | 2 |        | 1 |
        inc1 = temp; %  4 - 1        1 - 2        2 - 3
        temp = start2;
        start2 = start1;
        start1 = temp;
	end
end

%Output final values for loop increments to Attach Patch subroutine
if flag >= 1
    inc21 = minc+1;
    inc22 = 1;
    if swap == 1
%         temp = inc2;
%         inc2 = inc1;
%         inc1 = temp;
%         temp = start2;
%         start2 = start1;
%         start1 = temp;
        temp = inc22;
        inc22 = inc21;
        inc21 = temp;
    end
    start1 = start1 - 1;
    start2 = start2 - 1;
    inc21 = inc21*inc1;
    inc22 = inc22*inc2;
else
    return
end