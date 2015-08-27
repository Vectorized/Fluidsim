#ifndef Math_Lab_0_fastvec_h
#define Math_Lab_0_fastvec_h

#include <pmmintrin.h>
#include <math.h>
#include <iostream>

#define MI(m,index) ((float*)&m)[index]
#define MFLT(m) MI(m,0)

#define DEF_M_BASICS \
__m128 v;\
MV(){ v = _mm_set1_ps(0.f); }\
MV(float a) { v = _mm_set1_ps(a); }\
MV(const __m128 &v):v(v){};\
inline static float dot(const MV &a, const MV &b) {\
__m128 n0 = _mm_mul_ps(a.v, b.v);\
__m128 n1 = _mm_hadd_ps(n0, n0);\
__m128 r = _mm_hadd_ps(n1, n1);\
return MFLT(r);}\
inline float sum() {\
__m128 n0 = _mm_hadd_ps(v, v);\
__m128 r = _mm_hadd_ps(n0, n0);\
return MFLT(r);}\
inline void setMin(float f) {\
v = _mm_min_ps(v, _mm_set1_ps(f));}\
inline void setMax(float f) {\
v = _mm_max_ps(v, _mm_set1_ps(f));}\
inline static float distSq(const MV &a, const MV &b) { return dot(a, b); }\
inline float absSquared() { return dot(*this, *this); }\
inline float abs() const { return sqrtf(dot(*this, *this)); }\
inline float dot(const MV &o) { return dot(*this, o); }\
inline float operator[](int i) const { return ((float*)&v)[i]; }\
inline float & operator [](int i) {return ((float*)&v)[i]; }\
inline void normalize() { v = _mm_div_ps(v, _mm_set1_ps(abs())); }\
inline MV normalized() const { return MV(_mm_div_ps(v, _mm_set1_ps(abs()))); }\
inline float &x() { return MI(v, 0); }\
inline float &y() { return MI(v, 1); }\
inline float &z() { return MI(v, 2); }\
inline float &w() { return MI(v, 3); }\
inline float x() const { return MI(v, 0); }\
inline float y() const { return MI(v, 1); }\
inline float z() const { return MI(v, 2); }\
inline float w() const { return MI(v, 3); }\
inline void print();\
inline static MV cross(const MV &a, const MV &b) {\
__m128 r = _mm_sub_ps(\
_mm_mul_ps(a.v, _mm_shuffle_ps(b.v, b.v, _MM_SHUFFLE(3, 0, 2, 1))), \
_mm_mul_ps(b.v, _mm_shuffle_ps(a.v, a.v, _MM_SHUFFLE(3, 0, 2, 1)))); \
return MV(_mm_shuffle_ps(r, r, _MM_SHUFFLE(3, 0, 2, 1))); }\
inline operator const float*() const { return ((float*)&v); }\
inline operator float *() { return ((float*)&v); }\
inline MV &operator = (float f) { v = _mm_set1_ps(f); return *this; }\
inline MV floor() {\
__m128i t = _mm_cvtps_epi32(v);\
return MV(_mm_cvtepi32_ps(t));}\
void * operator new(size_t size) { return _mm_malloc(size, 1); }\
void operator delete(void *mem) { return _mm_free(mem); }

#define DEC_M_2 \
MV<2> xx(); MV<2> yy(); MV<2> zz();\
MV<2> xy(); MV<2> xz(); MV<2> yz();

#define DEF_M_2(dim) \
inline MV<2> MV<dim>::xx() { return MV<2>(MI(v, 0), MI(v, 0)); }\
inline MV<2> MV<dim>::yy() { return MV<2>(MI(v, 1), MI(v, 1)); }\
inline MV<2> MV<dim>::zz() { return MV<2>(MI(v, 2), MI(v, 2)); }\
inline MV<2> MV<dim>::xy() { return MV<2>(MI(v, 0), MI(v, 1)); }\
inline MV<2> MV<dim>::xz() { return MV<2>(MI(v, 0), MI(v, 2)); }\
inline MV<2> MV<dim>::yz() { return MV<2>(MI(v, 1), MI(v, 2)); }

#define DEC_M_3 \
MV<3> xyz(); MV<3> yzx(); MV<3> zxy();

#define DEF_M_3(dim) \
inline MV<3> MV<dim>::xyz() { return MV<3>(MI(v, 0), MI(v, 1), MI(v, 2)); }\
inline MV<3> MV<dim>::yzx() { return MV<3>(MI(v, 1), MI(v, 2), MI(v, 0)); }\
inline MV<3> MV<dim>::zxy() { return MV<3>(MI(v, 2), MI(v, 0), MI(v, 1)); }

template <int dim> struct MV { };

template <> struct MV<4> {
	MV(float a, float b, float c, float d=0.f) { v = _mm_set_ps(d,c,b,a); }
	MV(float *f) { v = _mm_set_ps(f[3],f[2],f[1],f[0]); }
	inline void set(float *f) { v = _mm_set_ps(f[3],f[2],f[1],f[0]); }
	DEF_M_BASICS;
	DEC_M_2;
	DEC_M_3;
};

template <> struct MV<3> {
	MV(float a, float b, float c) { v = _mm_set_ps(0.f,c,b,a); }
	MV(float *f) { v = _mm_set_ps(0.f,f[2],f[1],f[0]); }
	inline void set(float *f) { v = _mm_set_ps(0.f,f[2],f[1],f[0]); }
	DEF_M_BASICS;
	DEC_M_2;
	DEC_M_3;
};

template <> struct MV<2> {
	MV(float a, float b) { v = _mm_set_ps(0.f,0.f,b,a); }
	MV(float *f) { v = _mm_set_ps(0.f,0.f,f[1],f[0]); }
	inline void set(float *f) { v = _mm_set_ps(0.f,0.f,f[1],f[0]); }
	DEF_M_BASICS;
	DEC_M_2;
};

DEF_M_2(2); DEF_M_2(3); DEF_M_2(4);
DEF_M_3(3); DEF_M_3(4);

#define DEF_M_OP(op,op_inc,func)\
template <int dim> inline void operator op_inc (MV<dim> &b, const MV<dim> &o) { b.v = func(b.v, o.v); }\
template <int dim> inline MV<dim> operator op (const MV<dim> &b, const MV<dim> &o) { return MV<dim>(func(b.v, o.v)); }\
template <int dim> inline void operator op_inc (MV<dim> &b, float f) { b.v = func(b.v, _mm_set1_ps(f)); }\
template <int dim> inline MV<dim> operator op (const MV<dim> &b, float f) { return MV<dim>(func(b.v, _mm_set1_ps(f))); }\
template <int dim> inline MV<dim> operator op (float f, const MV<dim> &b) { return MV<dim>(func(b.v, _mm_set1_ps(f))); }
DEF_M_OP(+, +=, _mm_add_ps)
DEF_M_OP(-, -=, _mm_sub_ps)
DEF_M_OP(*, *=, _mm_mul_ps)
DEF_M_OP(/, /=, _mm_div_ps)

template <int dim> inline static MV<dim> lerp(const MV<dim>&a, const MV<dim> &b, float alpha) {
	return alpha * (a-b) + a;
}
template <int dim> inline static MV<dim> cubicInterpolate(const MV<dim> &p0, const MV<dim> &p1,
														  const MV<dim> &p2, const MV<dim> &p3, float t ) {
	MV<dim> p0p1 = lerp( p0, p1, t + 1 );
	MV<dim> p1p2 = lerp( p1, p2, t );
	MV<dim> p2p3 = lerp( p2, p3, t - 1 );
	MV<dim> p0p1_p1p2 = lerp( p0p1, p1p2, 0.5f * ( t + 1 ) );
	MV<dim> p1p2_p2p3 = lerp( p1p2, p2p3, 0.5f * t );
	return lerp( p0p1_p1p2, p1p2_p2p3, t );
}

template <int dim>inline bool operator == (const MV<dim> &a, const MV<dim> &b) {
	float r;
	_mm_store_ss(&r, _mm_cmpeq_ss(a.v, b.v));
	return r;
}
template <int dim> inline bool operator != (const MV<dim> &a, const MV<dim> &b) { return !(a == b); }

#define DEF_M_OSTREAM(dim,sprintf) \
inline ostream& operator<<(ostream& os, const MV<dim>& m) { char temp[64];sprintf;os << temp;return os; }
DEF_M_OSTREAM(4, sprintf(temp, "[ %.4f, %.4f, %.4f, %.4f ]",m[0],m[1],m[2],m[3]));
DEF_M_OSTREAM(3, sprintf(temp, "[ %.4f, %.4f, %.4f ]",m[0],m[1],m[2]));
DEF_M_OSTREAM(2, sprintf(temp, "[ %.4f, %.4f ]",m[0],m[1]));

#define DEF_PRINT(dim) inline void MV<dim>::print() { cout << *this << endl; }
DEF_PRINT(4);
DEF_PRINT(3);
DEF_PRINT(2);

template <int dim> inline MV<dim> max(const MV<dim> &a, const MV<dim> &b)
{
	return MV<dim>(_mm_max_ps(a.v, b.v));
}

template <int dim> inline MV<dim> min(const MV<dim> &a, const MV<dim> &b)
{
	return MV<dim>(_mm_min_ps(a.v, b.v));
}

template <int dim> inline float dot(const MV<dim> &a, const MV<dim> &b) { return MV<dim>::dot(a,b); }

template <int dim> inline MV<dim> cross(const MV<dim> &a, const MV<dim> &b) { return MV<dim>::cross(a,b); }

#include "mini_vec.h"
typedef MV<4> float4;
typedef mini_vec<MV<4>> float4v;

#endif
