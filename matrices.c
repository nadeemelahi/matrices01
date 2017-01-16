#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <assert.h>


// 180degrees = 3.14159265(PI)
#define PI 3.14159265
#define DegToRadxfactor PI/180
#define RadToDegxfactor 180/PI

#define degreesToRadians(x) (x*DegToRadxfactor)
#define radiansToDegrees(x) (x*RadToDegxfactor)

// row major 4x4 matrices 
// although opengl stores 
// in column major 
// 16 element matrix

typedef float M4[4][4];

void printM4(M4 m){
   printf("\n");
   //printf("size of row 0: %d\n", sizeof(m[0])/sizeof(float) ); // 4
   int i,j;
   printf("  ---                                      ---\n");
   for (i = 0; i<4; i++){
      printf("  |");
      for (j = 0; j<4; j++){
	 printf(" %8.3f ", m[i][j]);
      }
      printf("   |\n");
   }
   printf("  ---                                      ---\n");
   printf("\n");
}

void setM4ToIdentity(M4 m){
   int i,j;
   for (i = 0; i<4; i++){
      for (j = 0; j<4; j++){
	 m[i][j] = 0;
      }
   }
   m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.0;
}

void setM4ToTranslation(
      M4 m,
      float x, float y, float z ){

   setM4ToIdentity(m);

   //4th column 
   m[0][3] = x;
   m[1][3] = y;
   m[2][3] = z;
}

void setM4ToScale(
      M4 m,
      float sx, float sy, float sz){

   setM4ToIdentity(m);

   m[0][0] = sx;
   m[1][1] = sy;
   m[2][2] = sz;
}
void setM4ToRotateAboutX(
      M4 m,
      float degrees ){ 
   //Takes rotation is degrees not radians
   float radians = degreesToRadians(degrees);
   setM4ToIdentity(m);

   m[1][1] = cosf(radians);
   m[2][2] = m[1][1];

   m[1][2] = sinf(radians);
   m[2][1] = -m[1][2];

}

void setM4ToRotateAboutY(
      M4 m,
      float degrees ){ 
   //Takes rotation is degrees not radians
   float radians = degreesToRadians(degrees);
   setM4ToIdentity(m);

   m[0][0] = cosf(radians);
   m[2][2] = m[0][0];

   m[2][0] = sinf(radians);
   m[0][2] = -m[2][0];

}

void setM4ToRotateAboutZ(
      M4 m,
      float degrees ){ 
   //Takes rotation is degrees not radians
   float radians = degreesToRadians(degrees);
   setM4ToIdentity(m);

   m[0][0] = cosf(radians);
   m[1][1] = m[0][0];

   m[0][1] = sinf(radians);
   m[1][0] = -m[0][1];
}

void multiplyM4xM4(
      M4 a,
      M4 b,
      M4 c ){

   int rows, cols, step; 
   for (rows=0; rows<4; rows++){
      for (cols=0; cols<4; cols++){

	 c[rows][cols] = 0;

	 for(step=0; step<4; step++){
	    c[rows][cols] += a[rows][step]*b[step][cols];
	 }
      }
   }

}

//Euler rotation about arbitrary axis represented by unit vector and rotation amount(degrees/radians) Right Hand Rule
//http://www.cprogramming.com/tutorial/3d/quaternions.html
void setM4ToRotateAboutArbitraryVectorEuler(
      M4 m,
      float x, float y, float z, float degrees){

   setM4ToIdentity(m);

   float rads = degreesToRadians(degrees);

   float s = sinf(rads);
   float c = cosf(rads);
   float t = 1 - c;

   m[0][0] = t*x*x + c;
   m[0][1] = t*x*y + s*z;
   m[0][2] = t*x*z - s*y;
   
   m[1][0] = t*x*y - s*z;
   m[1][1] = t*y*y + c;
   m[1][2] = t*y*z + s*x;

   m[2][0] = t*x*z + s*y;
   m[2][1] = t*y*z - s*x;
   m[2][2] = t*z*z + c;

}

//http://physicsforgames.blogspot.ca/2010/02/quaternions.html
typedef struct Quaternion {
   float x;
   float y;
   float z;
   float w;
} Quatn;

void setQuatnToIdentity(Quatn q){
   q.x = q.y = q.z = 0.0;
   q.w = 1.0;
}
Quatn QuatFromAxisAngle(
      float x, float y, float z,
      float degrees){

     Quatn result;
     float rads = degreesToRadians(degrees);
     float s = sin(rads/2);
     float c = cos(rads/2);

     result.w = c;
     result.x = x * s;
     result.y = y * s;
     result.z = z * s;

     return result;
}

void setM4ToRotateAboutArbitraryVectorQuatn(
      M4 m,
      Quatn q ){

   setM4ToIdentity(m);

   //helper quantities - calculate these up front
   //to avoid redundancies
   float xSq = q.x * q.x;
   float ySq = q.y * q.y;
   float zSq = q.z * q.z;
   float wSq = q.w * q.w;
   float twoX = 2.0f * q.x;
   float twoY = 2.0f * q.y;
   float twoW = 2.0f * q.w;
   float xy = twoX * q.y;
   float xz = twoX * q.z;
   float yz = twoY * q.z;
   float wx = twoW * q.x;
   float wy = twoW * q.y;
   float wz = twoW * q.z;
   //fill in the first row
   m[0][0] = wSq + xSq - ySq - zSq;
   m[0][1] = xy + wz;
   m[0][2] = xz - wy;
   //fill in the second row
   m[1][0] = xy - wz;
   m[1][1] = wSq - xSq + ySq - zSq;
   m[1][2] = yz + wx;
   //fill in the third row
   m[2][0] = xz + wy;
   m[2][1] = yz - wx;
   m[2][2] = wSq - xSq - ySq + zSq;
}

Quatn QuatMultiply(Quatn q1, Quatn q2)
{
     Quatn result;
     result.w = q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z;
     result.x = q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y;
     result.y = q1.w*q2.y + q1.y*q2.w + q1.x*q2.z - q1.z*q2.x;
     result.z = q1.w*q2.z + q1.z*q2.w + q1.x*q2.y - q1.y*q2.x;
     return result;
}

Quatn QuatNormalize(Quatn q){
     Quatn result;
     float sq = q.x * q.x;
     sq += q.y * q.y;
     sq += q.z * q.z;
     sq += q.w * q.w;
     //detect badness
     //assert(sq > 0.1f);//#include <assert.h>
     float inv = 1.0f / sqrt(sq);
     result.x = q.x * inv;
     result.y = q.y * inv;
     result.z = q.z * inv;
     result.w = q.w * inv;
     return result;
}

//NLerp is a linear interpolation
Quatn QuatBlend(Quatn i, Quatn f, float blend){
     Quatn result;
     float dot = i.w*f.w + i.x*f.x + i.y*f.y + i.z*f.z;
     float blendI = 1.0f - blend;
     if(dot < 0.0f)
     {
          Quatn tmpF;
          tmpF.w = -f.w;
          tmpF.x = -f.x;
          tmpF.y = -f.y;
          tmpF.z = -f.z;
          result.w = blendI*i.w + blend*tmpF.w;
          result.x = blendI*i.x + blend*tmpF.x;
          result.y = blendI*i.y + blend*tmpF.y;
          result.z = blendI*i.z + blend*tmpF.z;
     }
     else
     {
          result.w = blendI*i.w + blend*f.w;
          result.x = blendI*i.x + blend*f.x;
          result.y = blendI*i.y + blend*f.y;
          result.z = blendI*i.z + blend*f.z;
     }
     result = QuatNormalize(result);
     return result;
}

// ortho projection matrix
//ex params (m4, 
//           0.001, 100.0, //near far
//           0.0, 320.0,   //left right
//           0.0, 480.0    //top bottom
//           )
void setM4ToFrustumOrtho( 
      M4 m,
      float near, float far, //create depth of field
      float left, float right, //rectangle view area
      float top, float bottom  ){

   m[0][0] = 2.0 / (right - left);
   m[1][0] = 0.0;
   m[2][0] = 0.0;
   m[3][0] = 0.0;

   m[0][1] = 0.0;
   m[1][1] = 2.0 / (top - bottom);
   m[2][1] = 0.0;
   m[3][1] = 0.0;
 
   m[0][2] = 0.0;
   m[1][2] = 0.0;
   m[2][2] = -2.0 / (far - near);
   m[3][2] = 0.0;
 
   m[0][3] = -(right + left) / (right - left);
   m[1][3] = -(top + bottom) / (top - bottom);
   m[2][3] = -(far + near) / (far - near);
   m[3][3] = 1.0;

}

// perspective projection matrix
// ex params(m4, 
//           0.001, 100.0, //near far
//           45.0, 0.75    //angleOfView aspectRatio
//           )
void setM4ToFrustumPerspective( 
      M4 m,
      float near, float far,
      float angleOfView,
      float aspectRatio ){

   float size = near * tanf(degreesToRadians(angleOfView));
   float left = -size;
   float right = size;
   float bottom = -size/aspectRatio;
   float top = size/aspectRatio;

   m[0][0] = 2 * near / (right - left);
   m[1][0] = 0.0;
   m[2][0] = 0.0;
   m[3][0] = 0.0;

   m[0][1] = 0.0;
   m[1][1] = 2 * near / (top - bottom);
   m[2][1] = 0.0;
   m[3][1] = 0.0;
 
   m[0][2] = (right + left) / (right - left);
   m[1][2] = (top + bottom) / (top - bottom);
   m[2][2] = -(far + near) / (far - near);
   m[3][2] = -1;
 
   m[0][3] = 0.0;
   m[1][3] = 0.0;
   m[2][3] = -(2 * far * near) / (far - near);
   m[3][3] = 0.0;


}

int main(){

   printf("********************************************\n");
   printf("******     hello M4[4][4]     **************\n");
   printf("* indentity, translation, scale, rotation **\n");
   printf("********************************************\n\n");

   M4 m4;
   setM4ToIdentity(m4);
   printf("****** identity ****************************\n");
   printM4(m4);
   
   M4 m4translation;
   setM4ToTranslation(m4translation, 10, 50, 2);
   printf("****** translation *************************\n");
   printM4(m4translation);

   M4 m4scale;
   setM4ToIdentity(m4scale);
   setM4ToScale(m4scale, 1.2, 3, 0.5);
   printf("****** scale *******************************\n");
   printM4(m4scale);
   
   M4 m4rotatex;
   setM4ToIdentity(m4rotatex);
   setM4ToRotateAboutX(m4rotatex, 60);
   printf("****** Rotate 60degrees about X  RHR********\n");
   printM4(m4rotatex);
   
   M4 m4rotatey;
   setM4ToIdentity(m4rotatey);
   setM4ToRotateAboutY(m4rotatey, 150);
   printf("****** Rotate 60degrees about Y  RHR********\n");
   printM4(m4rotatey);
   
   M4 m4rotatez;
   setM4ToIdentity(m4rotatez);
   setM4ToRotateAboutZ(m4rotatez, 10);
   printf("****** Rotate 60degrees about Z  RHR********\n");
   printM4(m4rotatez);


   printf("****** All above ***************************\n");
   printf("*** translation, scale, rotate xyz *********\n");
   printf("*** doing rotations 1 at a time ************\n");
   M4 m4scaletranslation;
   multiplyM4xM4(m4scale, m4translation, m4scaletranslation); //translation applied first
   M4 m4rotatexy;
   multiplyM4xM4(m4rotatex,m4rotatey,m4rotatexy);

   M4 m4rotatexyz;
   multiplyM4xM4(m4rotatez, m4rotatexy, m4rotatexyz);

   multiplyM4xM4(m4rotatexyz,m4scaletranslation,m4);
   
   printM4(m4);
  


   
   //begins matrix stuff
   printf("********************************************\n");
   printf("*****    hello Quaternions   ***************\n");
   printf("*** doing all 3 rotations at once **********\n");
   printf("*** using arbitrary vector and angle *******\n");
   printf("********************************************\n");
 
   printf("*** using Euler rotation equations *********\n");
   M4 m4RotateAboutVector;
   setM4ToIdentity(m4RotateAboutVector);
   setM4ToRotateAboutArbitraryVectorEuler(m4RotateAboutVector,0.2,1.1,0.5,25);
   printM4(m4RotateAboutVector);



   printf("*** using Quaternion rotation equations ****\n");
   Quatn q; 
   M4 mq;
   setQuatnToIdentity(q);

   q = QuatFromAxisAngle(0.2,1.1,0.5,25);
   printf("%f %f %f %f\n", q.w, q.x, q.y, q.z);
   
   setM4ToRotateAboutArbitraryVectorQuatn(mq,q);

   printM4(mq);
  
   exit(0);
}         
