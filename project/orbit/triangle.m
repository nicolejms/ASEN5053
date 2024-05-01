function [ang1 ang2 ang3 ] = triangle (s1,s2,s3)
%TRIANGLE calculates the degrees of all angles given the length of all sides.
% s1, s2, and s3 are the lengths of the sides of the triangle.
s= (s1+s2+s3)/ (2);
Area=sqrt (s* (s-s1)* (s-s2)* (s-s3));
if Area > 0
    ang1=acosd ( (s2^2 + s3^2 - s1^2)/ (2*s2*s3));
    ang2=acosd ( (s3^2 + s1^2 - s2^2)/ (2*s3*s1));
    ang3=acosd ( (s1^2 + s2^2 - s3^2)/ (2*s1*s2));
    fprintf ('The angle opposite side 1 is %.2f deg\\n',ang1);
    fprintf ('The angle opposite side 2 is %.2f deg\\n',ang2);
    fprintf ('The angle opposite side 3 is %.2f deg\\n',ang3);
elseif Area >= 0
    fprintf ('The 3 sides do NOT form an open triangle!\\n');
end
end
