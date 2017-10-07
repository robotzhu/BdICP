%由欧拉角到旋转矩阵的变换
% R = [r1 r2 r3; r4 r5 r6; r7 r8 r9];
% OulerAngle = [Rx, Ry Rz];
function [R] = OulerToRota(OulerAngle)
Rx = OulerAngle(1);    Ry = OulerAngle(2);     Rz = OulerAngle(3);
sx = sin(Rx);   cx = cos(Rx);
sy = sin(Ry);   cy = cos(Ry);
sz = sin(Rz);   cz = cos(Rz);
r1 = cy * cz; 
r2 = sx * sy * cz - cx * sz;
r3 = cx * sy * cz + sx * sz;
r4 = cy * sz;
r5 = sx * sy * sz + cx * cz;
r6 = cx * sy * sz - sx * cz;
r7 = -sy;
r8 = sx * cy;
r9 = cx * cy;
R = [r1 r2 r3
     r4 r5 r6
     r7 r8 r9];
 
      