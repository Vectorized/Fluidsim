#ifndef helpers
#define helpers

#ifdef _WIN32
#include "GL/freeglut.h"
#elif __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "vecmath.h"

// Just a bunch of macros and helper functions/classes I have written
// over the course of the past few 50.017 assignments.
// Lots of crazy metaprogramming here, pushing C++ to the limits.
// Many of the stuff can't be found online... =)
// e.g. the variable length MAX macro, or the VEC macro.

// START OF HELPERS -------------------------------------------------------------


#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880168872420969808
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#include <set>
#include <vector>
#include <sstream>
#include <iostream>
#include <map>
#include <algorithm>
#include <string>

using namespace std;

// Macro frenzy!!!
template<class T>struct __max{T m;__max<T>(T i){m=i;}__max&operator,(T i){if(m<i)m=i;return*this;}};
struct __maxmake{template<class T>__max<T>operator,(T i){return __max<T>(i);}};
#define MAX(...) ((__maxmake(),__VA_ARGS__).m)
template<class T>struct __min{T m;__min<T>(T i){m=i;}__min&operator,(T i){if(i<m)m=i;return*this;}};
struct __minmake{template<class T>__min<T>operator,(T i){return __min<T>(i);}};
#define MIN(...) ((__minmake(),__VA_ARGS__).m)
#define CLAMP(n,lower,upper) max(lower,min(n,upper))
#define SAME_SIGN(A,B) (((A)<0)==((B)<0))
#define SCREEN_REFRESH_DELAY 15

template<class T>inline T roundAway(T x){return x>0?ceil(x):floor(x);}

template<class T>inline void sort(vector<T>&vec){sort(vec.begin(),vec.end());};
template<class T>inline vector<T>sorted(vector<T>&vec){vector<T>v=vec;sort(v);return v;};

// Macro for for each in vs2010, since it does not support c++11 range based for each loops
// The multiple for loops is all just to make the variables within the scope of the foreach
// Dirty, but works.

#define FOREACH_WITH_INDEX(element,index,container) \
if(!container.empty())\
for(size_t __i=1,__j=1,index=0;__i&&__j;)\
for(auto &__c=container;__i;__i=0)\
for(auto element=*__c.begin();__j;__j=0)\
for(auto __it=__c.begin();__it!=__c.end();__it++,element=__it!=__c.end()?*__it:element,++index)

#define FOREACH(element,container) \
FOREACH_WITH_INDEX(element,__index,container)

#define FORINDEX(index,container) \
for(size_t index=0;index<container.size();++index)


template<class T>inline size_t __gfi(set<T>&s,long i){return 0;}
template<class T>inline size_t __gfi(vector<T>&v,long i){return i<0?(v.size()+i):MIN(i,v.size());}
template<class T>inline size_t __gli(set<T>&s,long i){return i<0?(s.size()+i+1):MIN(i+1,s.size());}
template<class T>inline size_t __gli(vector<T> &v,long i){return i<0?(v.size()+i+1):MIN(i+1,v.size());}
template<class T>inline typename set<T>::iterator __gi(set<T>&s,size_t i){return s.begin();}
template<class T>inline typename vector<T>::iterator __gi(vector<T>&v,size_t i){return v.begin()+i;}
// firstIndex and lastIndex form a closed range.
#define FOREACH_WITHIN_WITH_INDEX(element,firstIndex,lastIndex,index,container) \
if(!container.empty())\
for(bool __i=1,__j=1,__k=1;__i&&__j&&__k;)\
for(auto &__c=container;__i;__i=0)\
for(size_t index=__gfi(__c,firstIndex),\
__lastIndex=__gli(__c,lastIndex);\
__j&&index<__lastIndex;__j=0)\
for(auto element=*__gi(__c,index);__k;__k=0)\
for(auto __it=__gi(__c,index);\
index<__lastIndex;__it++,index++,element=index<__lastIndex?*__it:element)

#define FOREACH_WITHIN(element,firstIndex,lastIndex,container) \
FOREACH_WITHIN_WITH_INDEX(element,firstIndex,lastIndex,__index,container)

#define FOREACH_MAP(key,value,map) \
if(!map.empty())\
for(size_t __i=1,__j=1;__i&&__j;)\
for(auto key=map.begin()->first;__i;__i=0)\
for(auto &value=map.begin()->second;__j;__j=0)\
for(auto __it=map.begin();__it!=map.end();__it++,\
key=__it!=map.end()?__it->first:key,\
value=__it!=map.end()?__it->second:value)

// Fancy stuff for python style range loops...
// don't use in performance sensitive portions.
struct __rangeMaker{
int t;long b,e,s;
__rangeMaker(){t=-1;}
__rangeMaker&operator,(long i){++t;if(t==0){e=i;s=i>0?1:-1;}else if(t==1){b=e;e=i;s=e-b>0?1:-1;}
else if(t==2){s=i;if((b<e&&s>0)||(b>e&&s<0))e=((e-b)/s)*s+b;
else b=e=0;}return*this;}};

#define FORRANGE(index,...) \
for(bool __i=1;__i;)\
for(__rangeMaker __r=(__rangeMaker(),__VA_ARGS__);__i;__i=0)\
for(long index=__r.b;index!=__r.e;index+=__r.s)

typedef Vector2f V2;
typedef Vector3f V3;
typedef Vector4f V4;
typedef Matrix2f M2;
typedef Matrix3f M3;
typedef Matrix4f M4;
typedef Quat4f Q4;
typedef vector<Vector2f> V2v;
typedef vector<Vector3f> V3v;
typedef vector<Vector4f> V4v;

// Overloads << for vectors to become shorthand for push_back
template<class T,class V>vector<T>&operator<<(vector<T>&v,V e){v.push_back(e);return v;}
// Overloads << for sets to become shorthand for insert.
template<class T,class V>set<T>&operator<<(set<T> &s,V e){s.insert(e);return s;}

// Helper class and Macro to quick initialize vectors for vs2010 backwards compatibility.
template <class T> struct __vec{
	vector<T> v;
	__vec&operator,(T i){v.push_back(i);return*this;}
	__vec&operator,(vector<T>&_v){FOREACH(e,_v)v.push_back(e);return*this;}
	__vec&operator,(set<T>&_s){FOREACH(e, _s)v.push_back(e);return*this;}
};
struct __vecmake{template<class T>__vec<T>operator,(T i){return __vec<T>(),i;}};

#define VEC(...) ((__vecmake(),__VA_ARGS__).v)
#define VEC_INIT(name,...) auto name=VEC(__VA_ARGS__)
#define VEC_T(type,...) ((__vec<type>(),__VA_ARGS__).v)
#define VEC_T_INIT(name,type,...) auto name=VEC_T(type,__VA_ARGS__)

#define FOREACH_ARGS(element,...)\
for(size_t __i=1,__j=1;__i&&__j;)\
for(auto __c=VEC(__VA_ARGS__);__i;__i=0)\
for(auto element=*__c.begin();__j;__j=0)\
for(auto __it=__c.begin();__it!=__c.end();__it++,element=__it!=__c.end()?*__it:element)

// This macro is just pure magic...
#define FOREACH_ARGS_WITH_INDEX(element,index,...)\
for(size_t __i=1,__j=1,index=0;__i&&__j;)\
for(auto __c=VEC(__VA_ARGS__);__i;__i=0)\
for(auto element=*__c.begin();__j;__j=0)\
for(auto __it=__c.begin();__it!=__c.end();__it++,element=__it!=__c.end()?*__it:element,++index)

// Helper class and Macro to quick initialize sets for vs2010 backwards compaitbility
template <class T> struct __set{
	set<T> s;
	__set&operator,(T i){s.insert(i);return*this;}
	__set&operator,(set<T>&_s){FOREACH(e,_s)s.insert(e);return*this;}
	__set&operator,(vector<T> &_v){FOREACH(e,_v)s.insert(e);return*this;}
};
struct __setmake{template <class T> __set<T>operator,(T i){return __set<T>(),i;}};
#define SET(...) ((__setmake(),__VA_ARGS__).s)
#define SET_INIT(name,...) auto name=SET(__VA_ARGS__)
#define SET_T(type,...) ((__set<type>(),__VA_ARGS__).s)
#define SET_T_INIT(name,type,...) auto name=SET_T(type,__VA_ARGS__)

template <class T> vector<T> split(string s, string delims, T defaultValue) {
	vector<T> elems;
	stringstream ss;
	ss.str(s);
	string tokenString;
	stringstream sss;
	size_t prev = 0, pos;
	while ((pos = s.find_first_of(delims, prev)) != string::npos) {
		sss.str(s.substr(prev, pos-prev));
		T token;
		if (pos>prev && sss>>token) elems<<token;
		else elems<<defaultValue;
		prev = pos+1;
		sss.clear();
	}
	if (prev < s.length()) {
		T token;
		sss.str(s.substr(prev, string::npos));
		if (sss>>token)	elems<<token;
		else elems<<defaultValue;
	}
	return elems;
}

template <class T> vector<vector<T>> eachCons(vector<T> v, size_t chunkLength, size_t step) {
	vector<vector<T>> chunks;
	for (size_t i=0; i+chunkLength<=v.size(); i+=step) {
		vector<T> c;
		for (size_t j=0; j<chunkLength; ++j) {
			c.push_back(v[i+j]);
		}
		chunks.push_back(c);
	}
	return chunks;
}

#define VEC_WALK(statement,forVar,vec) \
for(size_t __i=0;__i<vec.size();++__i){auto forVar=vec[__i];vec[__i]=statement;}

struct __strjoiner{
	stringstream ss;
	template <class T> __strjoiner& operator,(T i) {
		ss<<i; return *this; }
};
#define STR(...) ((__strjoiner(),__VA_ARGS__).ss.str())

template <class T, class U>
inline bool contains(vector<T> const &v, U const &x) {
	return find(v.begin(), v.end(), x)!=v.end();
}
template <class T, class U>
inline bool contains(set<T> const &v, U const &x) {
	return v.find(x)!=v.end();
}
template <class T, class U, class V>
inline bool contains(map<T,V> const &v, U const &x) {
	return v.find(x)!=v.end();
}
#define CONTAINS(container,element) contains(container,element)

inline M4 differentiateBasis(M4 b) {
	return M4(b[4], b[8]*2,  b[12]*3, 0,
			  b[5], b[9]*2,  b[13]*3, 0,
			  b[6], b[10]*2, b[14]*3, 0,
			  b[7], b[11]*2, b[15]*3, 0);
}

#define __MATRIX4F_OP_DEF(overload,op) \
inline M4 operator overload (M4 a, M4 b) { \
return M4(a[0] op b[0], a[4] op b[4], a[8]  op b[8],  a[12] op b[12], \
a[1] op b[1], a[5] op b[5], a[9]  op b[9],  a[13] op b[13], \
a[2] op b[2], a[6] op b[6], a[10] op b[10], a[14] op b[14], \
a[3] op b[3], a[7] op b[7], a[11] op b[11], a[15] op b[15]);}\
inline M4 operator overload (M4 a, float b) { \
return M4(a[0] op b, a[4] op b, a[8]  op b, a[12] op b, \
a[1] op b, a[5] op b, a[9]  op b, a[13] op b, \
a[2] op b, a[6] op b, a[10] op b, a[14] op b, \
a[3] op b, a[7] op b, a[11] op b, a[15] op b);}

#define __MATRIX3F_OP_DEF(overload,op) \
inline M3 operator overload (M3 a, M3 b) { \
return M3(a[0] op b[0], a[3] op b[3], a[6] op b[6], \
a[1] op b[1], a[4] op b[4], a[7] op b[7], \
a[2] op b[2], a[5] op b[5], a[8] op b[8]);}\
inline M3 operator overload (M3 a, float b) { \
return M3(a[0] op b, a[3] op b, a[6] op b, \
a[1] op b, a[4] op b, a[7] op b, \
a[2] op b, a[5] op b, a[8] op b);}

#define __MATRIX2F_OP_DEF(overload,op) \
inline M2 operator overload (M2 a, M2 b) { \
return M2(a[0] op b[0], a[2] op b[2], \
a[1] op b[1], a[3] op b[3]);} \
inline M2 operator overload (M2 a, float b) { \
return M2(a[0] op b, a[2] op b, \
a[1] op b, a[3] op b);}

#define __MATRIX_OP_DEF(macro) \
macro(+,+);  macro(-,-);  macro(->*,*);  macro(/,/);

__MATRIX_OP_DEF(__MATRIX4F_OP_DEF);
__MATRIX_OP_DEF(__MATRIX3F_OP_DEF);
__MATRIX_OP_DEF(__MATRIX2F_OP_DEF);

inline M4 transposeSq(V4 k) {
	return M4(k[0]*k[0], k[0]*k[1], k[0]*k[2], k[0]*k[3],
			  k[1]*k[0], k[1]*k[1], k[1]*k[2], k[1]*k[3],
			  k[2]*k[0], k[2]*k[1], k[2]*k[2], k[2]*k[3],
			  k[3]*k[0], k[3]*k[1], k[3]*k[2], k[3]*k[3]);
}

inline M4 makeG(V4 v0, V4 v1, V4 v2, V4 v3) {
	return M4(v0, v1, v2, v3, true);
}
inline M4 makeG(V3 v0, V3 v1, V3 v2, V3 v3) {
	return makeG(V4(v0,0), V4(v1,0), V4(v2,0), V4(v3,0));
}
inline M4 makeG(V2 v0, V2 v1, V2 v2, V2 v3) {
	return makeG(V4(v0,0), V4(v1,0), V4(v2,0), V4(v3,0));
}

inline V3 planeNorm(V3 p0, V3 p1, V3 p2) {
	return V3::cross(p1-p0, p2-p0).normalized();}

inline V4 plane(V3 p0, V3 p1, V3 p2) {
	V3 n = planeNorm(p0, p1, p2);
	float z = V3::dot(p0, n);
	return V4(n[0], n[1], n[2], -z);
}

inline float angle(V3 a, V3 b, V3 c, bool degrees) {
	return acosf(V3::dot(a-b, b-c)) * (degrees ? 180/M_PI : 1);
}

inline void print() { printf("\n"); }
inline void print(M4 m) { m.print(); }
inline void print(M3 m) { m.print(); }
inline void print(M2 m) { m.print(); }
inline void print(V4 v) { v.print(); }
inline void print(V3 v) { v.print(); }
inline void print(V2 v) { v.print(); }
template <class T> inline void print(T v) { cout<<v<<endl; }

inline float distToLine(V3 a, V3 p0, V3 p1) {
	return V3::cross(a-p0, a-p1).abs()/(p1-p0).abs();
}


template <class T> void print(vector<T> &v) {
	printf("[ ");
	for (size_t i=0; i<v.size(); ++i) {
		cout << v[i];
		if (i<v.size()-1) printf(", ");
	}
	printf(" ]\n");
}

template <class T> void print(set<T> &v) {
	printf("[ ");
	for (size_t i=0; i<v.size(); ++i) {
		cout << v[i];
		if (i<v.size()-1) printf(", ");
	}
	printf(" ]\n");
}

struct Color3f {
	GLfloat r,g,b;
	Color3f(V4 _v) { r=_v[0];g=_v[1];b=_v[2]; }
	Color3f(V3 _v) { r=_v[0];g=_v[1];b=_v[2]; }
	Color3f(GLfloat _r, GLfloat _g, GLfloat _b) { r=_r;g=_g;b=_b; }
	Color3f(int hexValue) {
		r = ((hexValue >> 16) & 0xFF) / 255.0;
		g = ((hexValue >> 8) & 0xFF) / 255.0;
		b = ((hexValue) & 0xFF) / 255.0;
	}
	V3 toV3() { return V3(r,g,b); }
	V4 toV4(GLfloat a=1.0) { return V4(r,g,b,a); }
};

struct Color4f {
	GLfloat r,g,b,a;
	Color4f(GLfloat _r, GLfloat _g, GLfloat _b, GLfloat _a) { r=_r;g=_g;b=_b;a=_a; }
	Color4f(V4 _v) { r=_v[0];g=_v[1];b=_v[2]; }
	Color4f(V3 _v, GLfloat _a=1.0) { r=_v[0];g=_v[1];b=_v[2];a=_a; }
	Color4f(int hexValue) {
		r = ((hexValue >> 24) & 0xFF) / 255.0;
		g = ((hexValue >> 16) & 0xFF) / 255.0;
		b = ((hexValue >> 8) &  0xFF) / 255.0;
		a = ((hexValue) & 0xFF) / 255.0;
	}
	V4 toV4() { return V4(r,g,b,a); }
	V3 toV3() { return V3(r,g,b); }
};

inline void HSVtoRGB(float h, float s, float v, float &r, float &g, float &b){
	int i;
	float f,p,q,t;
	if(s==0) {r=g=b=v;return;}
	h /= 60;
	i = floor(h);
	f = h-i;
	p = v*(1-s);
	q = v*(1-s*f);
	t = v*(1-s*(1-f));
	switch(i) {
		case 0:r=v;g=t;b=p;break;
		case 1:r=q;g=v;b=p;break;
		case 2:r=p;g=v;b=t;break;
		case 3:r=p;g=q;b=v;break;
		case 4:r=t;g=p;b=v;break;
		default:r=v;g=p;b=q;break;
	}
}

struct TriNeighbors {
	size_t p0,p1;
	TriNeighbors(size_t _p0=0,size_t _p1=0){p0=_p0;p1=_p1;}
};

struct TriIndexes {
	size_t p0, p1, p2;
	TriIndexes(size_t _p0=0, size_t _p1=0, size_t _p2=0){p0=_p0;p1=_p1;p2=_p2;}
};

#define __APPROX_VEC_DEF(type)\
inline bool approx( const type &lhs, const type &rhs ) {\
const float eps = 1e-8f; \
return ( lhs - rhs ).absSquared() < eps;}

__APPROX_VEC_DEF(V4);
__APPROX_VEC_DEF(V3);
__APPROX_VEC_DEF(V2);

struct Face {
	size_t p[3];
	Face(size_t p0, size_t p1, size_t p2) {
		p[0] = p0; p[1] = p1; p[2] = p2;
	}
	Face() {}
	size_t& operator [] (size_t i) { return p[i]; };
};

inline V3 mid(V3 p0, V3 p1) {
	return (p0+p1)/2;
}

inline Q4 quatFrom2Vectors(V3 u, V3 v){
	float angle = acos(V3::dot(u.normalized(), v.normalized()));
	V3 w = V3::cross(u, v).normalized();
	Q4 q;
	q.setAxisAngle(angle, w);
	return q;
}

template <class T>
struct __typeconverter{
	stringstream ss;
	template <class V> T convert(V i) {
		T o;ss<<i;ss>>o;return o;
	}
};
#define CONVERT(type,var) (__typeconverter<type>().convert(var))
#define INT(var) CONVERT(int,var)
#define FLOAT(var) CONVERT(float,var)
#define LONG(var) CONVERT(long,var)

inline string getTillNextLine(istream &stream) {
	string s;char c;
	while (stream.peek()!='\n'){stream.get(c);s.push_back(c);}
	stream.get(c);return s;
}
inline string ltrim(string s,const char*t=" \t\n\r\f\v"){s.erase(0,s.find_first_not_of(t));return s;}
inline string rtrim(string s,const char*t=" \t\n\r\f\v"){s.erase(s.find_last_not_of(t)+1);return s;}
inline string trim(string s,const char*t=" \t\n\r\f\v"){return ltrim(rtrim(s,t),t);}

// Singleton for easy access to primitives here, there and everywhere.
class Globals {
private:
	map<string,bool>bools;
	map<string,long>longs;
	map<string,int>ints;
	map<string,float>floats;
	map<string,string>strings;
	map<string,double>doubles;
	Globals() {};
	Globals(Globals const&);
	void operator=(Globals const&);
	static Globals&i(){static Globals instance;return instance;}
public:
	static bool &getBool(string s){return i().bools[s];}
	static void setBool(string s, bool v){i().bools[s]=v;}
	static long &getLong(string s){return i().longs[s];}
	static void setLong(string s, long v){i().longs[s]=v;}
	static int  &getInt(string s){return i().ints[s];}
	static void setInt(string s, int v){i().ints[s]=v;}
	static float &getFloat(string s){return i().floats[s];}
	static void setFloat(string s, float v){i().floats[s]=v;}
	static string &getString(string s){return i().strings[s];}
	static void setString(string s, string v){i().strings[s]=v;}
	static double &getDouble(string s){return i().doubles[s];}
	static void setDouble(string s, double v){i().doubles[s]=v;}
};

inline bool startsWith(const string& haystack, const string& needle) {
	return needle.length() <= haystack.length()
	&& equal(needle.begin(), needle.end(), haystack.begin());
}

inline float randf() { return rand()/float(RAND_MAX); }

inline void calculateNorms(V3v &verts, V3v &norms, vector<Face> &faces) {
	norms.clear();
	norms.resize(verts.size());
	FOREACH(f, faces) {
		V3 n = planeNorm(verts[f[0]], verts[f[1]], verts[f[2]]);
		for (size_t i=0; i<3; ++i) norms[f[i]] += n;
	}
	for (size_t i=0; i<norms.size(); ++i) norms[i].normalize();
}

inline void drawFaces(V3v &verts, V3v &norms, vector<Face> &faces) {
	glBegin(GL_TRIANGLES);
	FOREACH(f, faces) {
		glNormal3fv(norms[f[0]]); glVertex3fv(verts[f[0]]);
		glNormal3fv(norms[f[1]]); glVertex3fv(verts[f[1]]);
		glNormal3fv(norms[f[2]]); glVertex3fv(verts[f[2]]);
	}
	glEnd();
}

inline string toLower(string &s) {
	string r = s;
	transform(r.begin(), r.end(), r.begin(), ::tolower);
	return r;
}

inline string toLower(char *c) {
	string s = c;
	return toLower(s);
}

inline string toUpper(string &s) {
	string r = s;
	transform(r.begin(), r.end(), r.begin(), ::toupper);
	return r;
}

inline string toUpper(char *c) {
	string s = c;
	return toUpper(s);
}

template <class T> inline T powi(T f, int pow) {
	T t = f;
	for (int i=1; i<pow; ++i) t*=f;
	return t;
}

inline int triIntersection(const V3 v0, const V3 v1, const V3 v2, const V3 origin,
						   const V3 direction, float &tOut, float epsilon = 0.000001)
{
	V3 e0 = v1 - v0;
	V3 e1 = v2 - v0;
	V3 P = V3::cross(direction, e1);
	float det = V3::dot(e0, P);
	if(det > -epsilon && det < epsilon) return 0;
	float invDet = 1.f / det;
	V3 T =  origin - v0;
 
	float u = V3::dot(T, P) * invDet;
	if(u < 0.f || u > 1.f) return 0;
 
	V3 Q = V3::cross(T, e0);
 
	float v = V3::dot(direction, Q) * invDet;
	if(v < 0.f || u + v  > 1.f) return 0;
 
	float t = V3::dot(e1, Q) * invDet;
 
	if(t > epsilon) {
		tOut = t;
		return 1;
	}
 
	return 0;
}

// END OF HELPERS -------------------------------------------------------------

#endif
