#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
struct Vec {        // Usage: time ./smallpt 5000 && xv image.ppm
  double x, y, z;                  // position, also color (r,g,b)
  Vec(double x_=0, double y_=0, double z_=0){ x=x_; y=y_; z=z_; }
  Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); }
  Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); }
  Vec operator*(double b) const { return Vec(x*b,y*b,z*b); }
  Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); }
  Vec& norm(){ return *this = *this * (1/sqrt(x*x+y*y+z*z)); }
  double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; } // cross:
  Vec operator%(Vec&b){return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
};
struct Ray { 
  Vec o, d; //origin, direction
  Ray(Vec o_, Vec d_) : o(o_), d(d_) {} //Makes o = first argument, etc (Constructor!)
};
struct RenderedRay { 
  Vec o, d, c; //origin, color
  RenderedRay(Vec o_=Vec(0,0,0), Vec d_=Vec(0,0,0), Vec c_=Vec(0,0,0)){ o=o_; d=d_; c=c_; }
  //RenderedRay(Vec o_, Vec d_, Vec c_) : o(o_), d(d_), c(c_) {} //Makes o = first argument, etc (Constructor!)
};
// o = origin, d = direction
// x = r.o + (r.d)*t (RAY INTERSECTION POINT!)
enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()
struct Sphere {
  double rad;       // radius
  Vec p, e, c;      // position, emission, color
  Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
  Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_):
    rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
  double intersect(const Ray &r) const { // returns distance, 0 if nohit
    Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    double t, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad;
    if (det<0) return 0; else det=sqrt(det);
    return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
  }
};
Sphere spheres[] = {//Scene: radius, position, emission, color, material
  Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
  Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght
  Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back
  Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(),           DIFF),//Frnt
  Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm
  Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
  Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirr
  Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glas
  Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), DIFF) //Lite
};
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }
inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); }
inline bool intersect(const Ray &r, double &t, int &id){
  double n=sizeof(spheres)/sizeof(Sphere), d, inf=t=1e20;
  for(int i=int(n);i--;) if((d=spheres[i].intersect(r))&&d<t){t=d;id=i;}
  return t<inf;
}
Vec radiance(const Ray &r, int depth, unsigned short *Xi){
	// Returning obj.e returns just the sphere's light emitting color
  double t;                               // distance to intersection
  int id=0;                               // id of intersected object
  if (!intersect(r, t, id)) return Vec(); // if miss, return black
  const Sphere &obj = spheres[id];        // obj: the hit object
  //obj.p, obj.e, obj.c = position, emission, color
  Vec x=r.o+r.d*t, n=(x-obj.p).norm(), nl=n.dot(r.d)<0?n:n*-1, f=obj.c; //obj.c is sphere color
  double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl
  // f = object color
  if (++depth>5) if (erand48(Xi)<p) f=f*(1/p); else return obj.e; //R.R.
  if (obj.refl == DIFF){                  // Ideal DIFFUSE reflection
    double r1=2*M_PI*erand48(Xi), r2=erand48(Xi), r2s=sqrt(r2);
    Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(), v=w%u;
    Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
    return obj.e + f.mult(radiance(Ray(x,d),depth,Xi));
  } else if (obj.refl == SPEC)            // Ideal SPECULAR reflection
    return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi));
  Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION
  bool into = n.dot(nl)>0;                // Ray from outside going in?
  double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t;
  if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
    return obj.e + f.mult(radiance(reflRay,depth,Xi));
  Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
  double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
  double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
  return obj.e + f.mult(depth>2 ? (erand48(Xi)<P ?   // Russian roulette
    radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
    radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
}
Vec sphereCam2(double delta, double columns, int multiplier) {
double x, y, z;
z = fmod(delta,columns); // z = height
double t = (delta - z); // t = 0-1 (Points around circle)
// sin(z*M_PI) = diameter of circle
z++;
x = cos((t*M_PI)/multiplier)*sin((z/columns));
y = sin((t*M_PI)/multiplier)*sin((z/columns));
//fprintf(stderr,"\r\nTest2: %f %f\r\n",sin((z/columns)*M_PI),z);
z = z - 1;
z = (z-columns/2)/(columns/2);
return Vec(x,y,z); // x, y, and z are from 0-1
}

Vec cubeCam(double delta, double columnz) {
double x, y, z;
double columns = columnz*4;
double il = floor(delta/columns); // i = row
double d = 2/columnz; // distance between columns
double l = delta - (il*columns); // leftover delta after vertical scan
double ol = floor(l/columnz);
double mult = delta - (ol*columnz) - (il*columns);
y = 1-(il*d);
if (l < (columnz)) {
x = (d*mult)-1;
z = 1;
} else if (l < (columnz*2)) {
x = 1;
z = 1-(d*mult);
} else if (l < (columnz*3)) {
x = 1-(d*mult);
z = -1;
} else if (l < (columnz*4)) {
x = -1;
z = (d*mult)-1;
}
//fprintf(stderr,"\r\n%f %f %f",mult,l,ol);
return Vec(x,y,z); // x, y, and z are from 0-1
}

Vec rayAngles(double delta, double columnz, Vec min, Vec max) {
double x, y, z;
double il = floor(delta/(columnz*columnz)); // i = row
double d = 2/columnz; // distance between columns
double l = delta - (il*columnz); // leftover delta after vertical scan
double ol = floor(l/columnz);
fprintf(stderr,"\r\n%f %f %f",il);
x = 1-(il*d);
y = 1-(il*d);
z = 1-(il*d);
double mult = delta - (ol*columnz) - (il*columnz);
return Vec(x,y,z);
}
int mainold(int argc, char *argv[]){
  int w=1024, h=768, samps = argc==2 ? atoi(argv[1])/4 : 1; // # samples
  Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir
  Vec cx=Vec(w*.5135/h), cy=(cx%cam.d).norm()*.5135, r, *c=new Vec[w*h];
//#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
  for (int y=0; y<h; y++){                       // Loop over image rows
    fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
    for (unsigned short x=0, Xi[3]={0,0,y*y*y}; x<w; x++)   // Loop cols
      for (int sy=0, i=(h-y-1)*w+x; sy<2; sy++)     // 2x2 subpixel rows
        for (int sx=0; sx<2; sx++, r=Vec()){        // 2x2 subpixel cols
          for (int s=0; s<samps; s++){
            double r1=2*erand48(Xi), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
            double r2=2*erand48(Xi), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);
            Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
                    cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;
            r = r + radiance(Ray(cam.o+d*140,d.norm()),0,Xi)*(1./samps);
          } // Camera rays are pushed ^^^^^ forward to start in interior
          c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25;
        }
  }
  FILE *f = fopen("image.ppm", "w");         // Write image to PPM file.
  fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
  for (int i=0; i<w*h; i++)
    fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}
bool vecEquals(Vec a, Vec b) {
return ((a.x == b.x) && (a.y == b.y) && (a.z == b.z));
}
int mainoldold(int argc, char *argv[]){
	int columns = argc==2 ? atoi(argv[1]) : 5;
	int multiplier = columns*columns;
	int sq = multiplier * multiplier;
	RenderedRay *c = new RenderedRay[multiplier*multiplier];
	for (int i=0; i < multiplier; i++) {
		unsigned short Xi[3]={0,0,i*i*i};
		Vec st = cubeCam(i,columns);
		for (int p=0; p < multiplier; p++) {
			Vec sp = cubeCam(p,columns);
			c[(i*multiplier)+p] = RenderedRay(st,sp,radiance(Ray(Vec(st.x*200,st.y*200,st.z*200) + Vec(50,52,295.6),sp),0,Xi));
		}
		fprintf(stderr,"\rRendering (Multiplier: %d) %5.2f%%",multiplier,100.*i/(multiplier - 1));
	}
	FILE *f = fopen("sphere.prm", "w"); // Write image to Portable Raymap file.
	fprintf(f, "R1\n%d", multiplier); // Line 1: Portable Raymap Type Identifier, Line 2: Portable Raymap Ray Multiplier
	Vec last = Vec(9,9,9);
	for (int i=0; i<multiplier*multiplier; i++){
			if (vecEquals(c[i].o,last) == false){
				last = c[i].o;
				fprintf(f,"\n%f %f %f\n", c[i].o.x, c[i].o.y, c[i].o.z/columns);
			}
		fprintf(f,"%f %f %f %d %d %d ", c[i].d.x, c[i].d.y, c[i].d.z, toInt(c[i].c.x), toInt(c[i].c.y), toInt(c[i].c.z));
	}
}
int main(int argc, char *argv[]) {
// Each face of cube has columns^2 points on it)
int columns = argc==2 ? atoi(argv[1]) : 5;
	int multiplier = columns*columns*4;
	int max = columns*columns*6;
	Vec *points = new Vec[columns*columns*6];
for (int i=0; i < columns*columns; i++) {
	unsigned short Xi[3]={0,0,i*i*i};
	Vec st = rayAngles(i,columns,Vec(-1,-1,-1),Vec(1,1,1));
	fprintf(stderr,"\r\nRendering (Multiplier: %d) %f %f %f",i,st.x,st.y,st.z);
	}
}
