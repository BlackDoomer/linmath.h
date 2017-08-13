#ifndef LINMATH_H
#define LINMATH_H

#include <math.h>
#include <float.h>

#ifdef _MSC_VER
	#define inline __inline
#endif

////////////////////////////////////////////////////////////////////////////////
// Implementation macros
////////////////////////////////////////////////////////////////////////////////

#ifndef LINMATH_H_NEAR_ZERO
	#define LINMATH_H_NEAR_ZERO FLT_EPSILON
#endif

/*
  Next macros are provided to implement optional row-major mode for matrices:
  __LH_M4E - access mat4x4 element
  __LH_POS - positive if not row-major, otherwise negative
  __LH_NEG - same as __LH_POS, but vice versa
*/

#ifndef LINMATH_H_ROW_MAJOR
	#define __LH_M4E(m,c,r) (m[c][r])
	#define __LH_POS(x) (x)
	#define __LH_NEG(x) (-x)
#else
	#define __LH_M4E(m,c,r) (m[r][c])
	#define __LH_POS(x) (-x)
	#define __LH_NEG(x) (x)
#endif

/*
  Union vec#_t types support unnamed vector elements access when
  supported by compiler - that is, for example,
    vec#_t.x = 1.0f;
  instead of
    vec#_t.i.x = 1.0f;
  Note that the named syntax (with ".i") is provided in any way for those who
  don't want to use less portable unnamed syntax.
*/

/* Unnamed structs - since C11; as compiler extensions in older standards. */
#if !defined(LINMATH_H_NO_UNNAMED_STRUCTS) && defined(__STDC_VERSION__)
	#if (__STDC_VERSION__ < 201112L) && \
	  ( !defined(__GNUC__) || defined(__STRICT_ANSI__) )
		#define LINMATH_H_NO_UNNAMED_STRUCTS
	#endif
#endif

#ifndef LINMATH_H_NO_UNNAMED_STRUCTS
	#define __LH_VECTOR_STRUCT(...) \
		struct {float __VA_ARGS__;}; \
		struct {float __VA_ARGS__;} i
#else
	#define __LH_VECTOR_STRUCT(...) \
		struct {float __VA_ARGS__;} i
#endif

#define LINMATH_H_DEFINE_VEC(n, ...) \
	typedef float vec##n[n]; \
	typedef union { \
		vec##n v; \
		__LH_VECTOR_STRUCT(__VA_ARGS__); \
	} vec##n##_t; \
	static inline void vec##n##_dup(vec##n r, vec##n v) \
	{ \
		for(int i=0; i<n; ++i) \
			r[i] = v[i]; \
	} \
	static inline void vec##n##_add(vec##n r, vec##n a, vec##n b) \
	{ \
		for(int i=0; i<n; ++i) \
			r[i] = a[i] + b[i]; \
	} \
	static inline void vec##n##_sub(vec##n r, vec##n a, vec##n b) \
	{ \
		for(int i=0; i<n; ++i) \
			r[i] = a[i] - b[i]; \
	} \
	static inline void vec##n##_scale(vec##n r, vec##n v, float s) \
	{ \
		for(int i=0; i<n; ++i) \
			r[i] = v[i] * s; \
	} \
	static inline float vec##n##_mul_inner(vec##n a, vec##n b) \
	{ \
		float p = 0.f; \
		for(int i=0; i<n; ++i) \
			p += b[i] * a[i]; \
		return p; \
	} \
	static inline float vec##n##_len(vec##n v) \
	{ \
		return sqrtf(vec##n##_mul_inner(v,v)); \
	} \
	static inline void vec##n##_norm(vec##n r, vec##n v) \
	{ \
		float k = 1.f / vec##n##_len(v); \
		vec##n##_scale(r, v, k); \
	} \
	static inline void vec##n##_reflect(vec##n r, vec##n v, vec##n o) \
	{ \
		float p = 2.f * vec##n##_mul_inner(v, o); \
		for(int i=0; i<n; ++i) \
			r[i] = v[i] - p*o[i]; \
	} \
	static inline void vec##n##_min(vec##n r, vec##n a, vec##n b) \
	{ \
		for(int i=0; i<n; ++i) \
			r[i] = (a[i] < b[i]) ? a[i] : b[i]; \
	} \
	static inline void vec##n##_max(vec##n r, vec##n a, vec##n b) \
	{ \
		for(int i=0; i<n; ++i) \
			r[i] = (a[i] > b[i]) ? a[i] : b[i]; \
	} \
	extern void __THIS_REQUIRES_A_SEMICOLON(void)
// LINMATH_H_DEFINE_VEC

////////////////////////////////////////////////////////////////////////////////
// Vectors types: vec2 / vec2_t, vec3 / vec3_t, vec4 / vec4_t
////////////////////////////////////////////////////////////////////////////////

LINMATH_H_DEFINE_VEC(2, x, y);
LINMATH_H_DEFINE_VEC(3, x, y, z);
LINMATH_H_DEFINE_VEC(4, x, y, z, w);

#ifndef LINMATH_H_NO_BUILDERS
	/*
	  NOTE: There's no macros for plain vectors constant initializers.
	  Use default C syntax for this purpose:
	    vec2 v_xy = {x,y};
	    vec3 v_xyz = {x,y,z};
	    vec4 v_xyzw = {x,y,z,w};
	*/

	/*
	  Function arguments builders for plain vectors.
	  Usage example:
	    void f(vec2 v_xy);
	    f(VEC2(1.0f, 2.0f));
	*/
	#define VEC2(x,y)     (vec2){x,y}
	#define VEC3(x,y,z)   (vec3){x,y,z}
	#define VEC4(x,y,z,w) (vec4){x,y,z,w}

	/*
	  Union types constant initializers, for compile-time initialization.
	  Usage example:
	    vec2_t vt_xy = iVEC2T(1.0f, 2.0f);
	*/
	#define iVEC2T(x,y)     {.v={x,y}}
	#define iVEC3T(x,y,z)   {.v={x,y,z}}
	#define iVEC4T(x,y,z,w) {.v={x,y,z,w}}

	/*
	  Function arguments builders for union types.
	  Usage example:
	    void f(vec2_t vt_xy);
	    f(VEC2T(1.0f, 2.0f));
	*/
	#define VEC2T(x,y)     (vec2_t)iVEC2T(x,y)
	#define VEC3T(x,y,z)   (vec3_t)iVEC3T(x,y,z)
	#define VEC4T(x,y,z,w) (vec4_t)iVEC4T(x,y,z,w)

	/*
	  Macros for casting plain vectors to union types.
	  Usage example:
	    vec2 v_xy = {1.0f, 2.0f};
	    vec2_t vt_xy = cVEC2T(v_xy);
	*/
	#define cVEC2T(v2) VEC2T(v2[0],v2[1])
	#define cVEC3T(v3) VEC3T(v3[0],v3[1],v3[2])
	#define cVEC4T(v4) VEC4T(v4[0],v4[1],v4[2],v4[3])
#endif

static inline void vec3_mul_cross(vec3 r, vec3 a, vec3 b)
{
	r[0] = a[1]*b[2] - a[2]*b[1];
	r[1] = a[2]*b[0] - a[0]*b[2];
	r[2] = a[0]*b[1] - a[1]*b[0];
}

static inline void vec4_mul_cross(vec4 r, vec4 a, vec4 b)
{
	vec3_mul_cross(r,a,b);
	r[3] = 1.f;
}

////////////////////////////////////////////////////////////////////////////////
// 4x4 matrices type: mat4x4
////////////////////////////////////////////////////////////////////////////////

typedef vec4 mat4x4[4];

static inline void mat4x4_identity(mat4x4 M)
{
	for(int i=0; i<4; ++i)
		for(int j=0; j<4; ++j)
			M[i][j] = (i==j) ? 1.f : 0.f;
}
static inline void mat4x4_dup(mat4x4 M, mat4x4 N)
{
	for(int i=0; i<4; ++i)
		for(int j=0; j<4; ++j)
			M[i][j] = N[i][j];
}
static inline void mat4x4_row(vec4 r, mat4x4 M, int i)
{
	for(int k=0; k<4; ++k)
		r[k] = __LH_M4E(M,k,i);
}
static inline void mat4x4_col(vec4 r, mat4x4 M, int i)
{
	for(int k=0; k<4; ++k)
		r[k] = __LH_M4E(M,i,k);
}
static inline void mat4x4_set_row(mat4x4 M, vec4 v, int i)
{
	for(int k=0; k<4; ++k)
		__LH_M4E(M,k,i) = v[k];
}
static inline void mat4x4_set_col(mat4x4 M, vec4 v, int i)
{
	for(int k=0; k<4; ++k)
		__LH_M4E(M,i,k) = v[k];
}
static inline void mat4x4_transpose(mat4x4 M, mat4x4 N)
{
	for(int j=0; j<4; ++j)
		for(int i=0; i<4; ++i)
			M[i][j] = N[j][i];
}
static inline void mat4x4_add(mat4x4 M, mat4x4 a, mat4x4 b)
{
	for(int i=0; i<4; ++i)
		vec4_add(M[i], a[i], b[i]);
}
static inline void mat4x4_sub(mat4x4 M, mat4x4 a, mat4x4 b)
{
	for(int i=0; i<4; ++i)
		vec4_sub(M[i], a[i], b[i]);
}
static inline void mat4x4_scale(mat4x4 M, mat4x4 a, float k)
{
	for(int i=0; i<4; ++i)
		vec4_scale(M[i], a[i], k);
}
static inline void mat4x4_scale_aniso(mat4x4 M, mat4x4 a, float x, float y, float z)
{
	vec4 a0, a1, a2;
	mat4x4_col(a0, a, 0);
	mat4x4_col(a1, a, 1);
	mat4x4_col(a2, a, 2);

	vec4_scale(a0, a0, x);
	vec4_scale(a1, a1, y);
	vec4_scale(a2, a2, z);

	mat4x4_set_col(M, a0, 0);
	mat4x4_set_col(M, a1, 1);
	mat4x4_set_col(M, a2, 2);

	for(int i=0; i<4; ++i)
		__LH_M4E(M,3,i) = __LH_M4E(a,3,i);
}
static inline void mat4x4_mul(mat4x4 M, mat4x4 a, mat4x4 b)
{
	mat4x4 temp;
	for(int c=0; c<4; ++c) for(int r=0; r<4; ++r) {
		__LH_M4E(temp,c,r) = 0.f;
		for(int k=0; k<4; ++k)
			__LH_M4E(temp,c,r) += __LH_M4E(a,k,r) * __LH_M4E(b,c,k);
	}
	mat4x4_dup(M, temp);
}
static inline void mat4x4_mul_vec4(vec4 r, mat4x4 M, vec4 v)
{
	for(int j=0; j<4; ++j) {
		r[j] = 0.f;
		for(int i=0; i<4; ++i)
			r[j] += __LH_M4E(M,i,j) * v[i];
	}
}
static inline void mat4x4_translate(mat4x4 T, float x, float y, float z)
{
	mat4x4_identity(T);
	__LH_M4E(T,3,0) = x;
	__LH_M4E(T,3,1) = y;
	__LH_M4E(T,3,2) = z;
}
static inline void mat4x4_translate_in_place(mat4x4 M, float x, float y, float z)
{
	vec4 t = {x, y, z, 0};
	vec4 r;
	for(int i=0; i<4; ++i) {
		mat4x4_row(r, M, i);
		__LH_M4E(M,3,i) += vec4_mul_inner(r, t);
	}
}
static inline void mat4x4_from_vec3_mul_outer(mat4x4 M, vec3 a, vec3 b)
{
	for(int i=0; i<4; ++i)
		for(int j=0; j<4; ++j)
			__LH_M4E(M,i,j) = ((i<3) && (j<3)) ? (a[i] * b[j]) : 0.f;
}
static inline void mat4x4_rotate(mat4x4 R, mat4x4 M, float x, float y, float z, float angle)
{
	float s = sinf(angle);
	float c = cosf(angle);
	vec3 u = {x, y, z};

	if(vec3_len(u) > LINMATH_H_NEAR_ZERO) {
		vec3_norm(u, u);
		mat4x4 T;
		mat4x4_from_vec3_mul_outer(T, u, u);

		mat4x4 S = {
			{            0.f, __LH_POS(u[2]), __LH_NEG(u[1]), 0.f },
			{ __LH_NEG(u[2]),            0.f, __LH_POS(u[0]), 0.f },
			{ __LH_POS(u[1]), __LH_NEG(u[0]),            0.f, 0.f },
			{            0.f,            0.f,            0.f, 0.f }
		};

		mat4x4_scale(S, S, s);

		mat4x4 C;
		mat4x4_identity(C);
		mat4x4_sub(C, C, T);

		mat4x4_scale(C, C, c);

		mat4x4_add(T, T, C);
		mat4x4_add(T, T, S);

		T[3][3] = 1.f;
		mat4x4_mul(R, M, T);
	} else {
		mat4x4_dup(R, M);
	}
}
static inline void mat4x4_rotate_X(mat4x4 Q, mat4x4 M, float angle)
{
	float s = sinf(angle);
	float c = cosf(angle);
	mat4x4 R = {
		{ 1.f,         0.f,         0.f, 0.f },
		{ 0.f,           c, __LH_POS(s), 0.f },
		{ 0.f, __LH_NEG(s),           c, 0.f },
		{ 0.f,         0.f,         0.f, 1.f }
	};
	mat4x4_mul(Q, M, R);
}
static inline void mat4x4_rotate_Y(mat4x4 Q, mat4x4 M, float angle)
{
	float s = sinf(angle);
	float c = cosf(angle);
	mat4x4 R = {
		{           c, 0.f, __LH_POS(s), 0.f },
		{         0.f, 1.f,         0.f, 0.f },
		{ __LH_NEG(s), 0.f,           c, 0.f },
		{         0.f, 0.f,         0.f, 1.f }
	};
	mat4x4_mul(Q, M, R);
}
static inline void mat4x4_rotate_Z(mat4x4 Q, mat4x4 M, float angle)
{
	float s = sinf(angle);
	float c = cosf(angle);
	mat4x4 R = {
		{           c, __LH_POS(s), 0.f, 0.f },
		{ __LH_NEG(s),           c, 0.f, 0.f },
		{         0.f,         0.f, 1.f, 0.f },
		{         0.f,         0.f, 0.f, 1.f }
	};
	mat4x4_mul(Q, M, R);
}
static inline void mat4x4_invert(mat4x4 T, mat4x4 M)
{
	float s[6];
	float c[6];
	s[0] = M[0][0]*M[1][1] - M[1][0]*M[0][1];
	s[1] = M[0][0]*M[1][2] - M[1][0]*M[0][2];
	s[2] = M[0][0]*M[1][3] - M[1][0]*M[0][3];
	s[3] = M[0][1]*M[1][2] - M[1][1]*M[0][2];
	s[4] = M[0][1]*M[1][3] - M[1][1]*M[0][3];
	s[5] = M[0][2]*M[1][3] - M[1][2]*M[0][3];

	c[0] = M[2][0]*M[3][1] - M[3][0]*M[2][1];
	c[1] = M[2][0]*M[3][2] - M[3][0]*M[2][2];
	c[2] = M[2][0]*M[3][3] - M[3][0]*M[2][3];
	c[3] = M[2][1]*M[3][2] - M[3][1]*M[2][2];
	c[4] = M[2][1]*M[3][3] - M[3][1]*M[2][3];
	c[5] = M[2][2]*M[3][3] - M[3][2]*M[2][3];

	/* Assumes it is invertible */
	float idet = 1.f/( s[0]*c[5]-s[1]*c[4]+s[2]*c[3]+s[3]*c[2]-s[4]*c[1]+s[5]*c[0] );

	T[0][0] = ( M[1][1] * c[5] - M[1][2] * c[4] + M[1][3] * c[3]) * idet;
	T[0][1] = (-M[0][1] * c[5] + M[0][2] * c[4] - M[0][3] * c[3]) * idet;
	T[0][2] = ( M[3][1] * s[5] - M[3][2] * s[4] + M[3][3] * s[3]) * idet;
	T[0][3] = (-M[2][1] * s[5] + M[2][2] * s[4] - M[2][3] * s[3]) * idet;

	T[1][0] = (-M[1][0] * c[5] + M[1][2] * c[2] - M[1][3] * c[1]) * idet;
	T[1][1] = ( M[0][0] * c[5] - M[0][2] * c[2] + M[0][3] * c[1]) * idet;
	T[1][2] = (-M[3][0] * s[5] + M[3][2] * s[2] - M[3][3] * s[1]) * idet;
	T[1][3] = ( M[2][0] * s[5] - M[2][2] * s[2] + M[2][3] * s[1]) * idet;

	T[2][0] = ( M[1][0] * c[4] - M[1][1] * c[2] + M[1][3] * c[0]) * idet;
	T[2][1] = (-M[0][0] * c[4] + M[0][1] * c[2] - M[0][3] * c[0]) * idet;
	T[2][2] = ( M[3][0] * s[4] - M[3][1] * s[2] + M[3][3] * s[0]) * idet;
	T[2][3] = (-M[2][0] * s[4] + M[2][1] * s[2] - M[2][3] * s[0]) * idet;

	T[3][0] = (-M[1][0] * c[3] + M[1][1] * c[1] - M[1][2] * c[0]) * idet;
	T[3][1] = ( M[0][0] * c[3] - M[0][1] * c[1] + M[0][2] * c[0]) * idet;
	T[3][2] = (-M[3][0] * s[3] + M[3][1] * s[1] - M[3][2] * s[0]) * idet;
	T[3][3] = ( M[2][0] * s[3] - M[2][1] * s[1] + M[2][2] * s[0]) * idet;
}
static inline void mat4x4_orthonormalize(mat4x4 R, mat4x4 M)
{
	float s = 1.f;
	vec3 h;

	vec4 r0, r1, r2, _r3;
	mat4x4_col(r0,  M, 0);
	mat4x4_col(r1,  M, 1);
	mat4x4_col(r2,  M, 2);
	mat4x4_col(_r3, M, 3);

	vec3_norm(r2, r2);

	s = vec3_mul_inner(r1, r2);
	vec3_scale(h, r2, s);
	vec3_sub(r1, r1, h);
	vec3_norm(r2, r2);

	s = vec3_mul_inner(r1, r2);
	vec3_scale(h, r2, s);
	vec3_sub(r1, r1, h);
	vec3_norm(r1, r1);

	s = vec3_mul_inner(r0, r1);
	vec3_scale(h, r1, s);
	vec3_sub(r0, r0, h);
	vec3_norm(r0, r0);

	mat4x4_set_col(R, r0,  0);
	mat4x4_set_col(R, r1,  1);
	mat4x4_set_col(R, r2,  2);
	mat4x4_set_col(R, _r3, 3);
}

static inline void mat4x4_frustum(mat4x4 M, float l, float r, float b, float t, float n, float f)
{
	__LH_M4E(M,0,0) = 2.f*n / (r-l);
	__LH_M4E(M,0,1) = __LH_M4E(M,0,2) = __LH_M4E(M,0,3) = 0.f;

	__LH_M4E(M,1,1) = 2.f*n / (t-b);
	__LH_M4E(M,1,0) = __LH_M4E(M,1,2) = __LH_M4E(M,1,3) = 0.f;

	__LH_M4E(M,2,0) = (r+l) / (r-l);
	__LH_M4E(M,2,1) = (t+b) / (t-b);
	__LH_M4E(M,2,2) = -(f+n) / (f-n);
	__LH_M4E(M,2,3) = -1.f;

	__LH_M4E(M,3,2) = -2.f * (f*n) / (f-n);
	__LH_M4E(M,3,0) = __LH_M4E(M,3,1) = __LH_M4E(M,3,3) = 0.f;
}
static inline void mat4x4_ortho(mat4x4 M, float l, float r, float b, float t, float n, float f)
{
	__LH_M4E(M,0,0) = 2.f / (r-l);
	__LH_M4E(M,0,1) = __LH_M4E(M,0,2) = __LH_M4E(M,0,3) = 0.f;

	__LH_M4E(M,1,1) = 2.f / (t-b);
	__LH_M4E(M,1,0) = __LH_M4E(M,1,2) = __LH_M4E(M,1,3) = 0.f;

	__LH_M4E(M,2,2) = -2.f / (f-n);
	__LH_M4E(M,2,0) = __LH_M4E(M,2,1) = __LH_M4E(M,2,3) = 0.f;

	__LH_M4E(M,3,0) = -(r+l) / (r-l);
	__LH_M4E(M,3,1) = -(t+b) / (t-b);
	__LH_M4E(M,3,2) = -(f+n) / (f-n);
	__LH_M4E(M,3,3) = 1.f;
}
static inline void mat4x4_perspective(mat4x4 m, float y_fov, float aspect, float n, float f)
{
	/*
	  NOTE: Degrees are an unhandy unit to work with.
	  linmath.h uses radians for everything!
	*/

	float const a = 1.f / tanf(y_fov/2.f);

	__LH_M4E(m,0,0) = a / aspect;
	__LH_M4E(m,0,1) = 0.f;
	__LH_M4E(m,0,2) = 0.f;
	__LH_M4E(m,0,3) = 0.f;

	__LH_M4E(m,1,0) = 0.f;
	__LH_M4E(m,1,1) = a;
	__LH_M4E(m,1,2) = 0.f;
	__LH_M4E(m,1,3) = 0.f;

	__LH_M4E(m,2,0) = 0.f;
	__LH_M4E(m,2,1) = 0.f;
	__LH_M4E(m,2,2) = -((f+n)/(f-n));
	__LH_M4E(m,2,3) = -1.f;

	__LH_M4E(m,3,0) = 0.f;
	__LH_M4E(m,3,1) = 0.f;
	__LH_M4E(m,3,2) = -((2.f*f*n)/(f-n));
	__LH_M4E(m,3,3) = 0.f;
}
static inline void mat4x4_look_at(mat4x4 m, vec3 eye, vec3 center, vec3 up)
{
	/*
	  Adapted from Android's OpenGL Matrix.java.
	  See the OpenGL GLUT documentation for gluLookAt for a description
	  of the algorithm. We implement it in a straightforward way.

	  TODO: The negation of vector can be spared by swapping the order of
	  operands in the following cross products in the right way.
	*/

	vec3 f;
	vec3_sub(f, center, eye);
	vec3_norm(f, f);

	vec3 s;
	vec3_mul_cross(s, f, up);
	vec3_norm(s, s);

	vec3 t;
	vec3_mul_cross(t, s, f);

	__LH_M4E(m,0,0) =  s[0];
	__LH_M4E(m,0,1) =  t[0];
	__LH_M4E(m,0,2) = -f[0];
	__LH_M4E(m,0,3) =  0.f;

	__LH_M4E(m,1,0) =  s[1];
	__LH_M4E(m,1,1) =  t[1];
	__LH_M4E(m,1,2) = -f[1];
	__LH_M4E(m,1,3) =  0.f;

	__LH_M4E(m,2,0) =  s[2];
	__LH_M4E(m,2,1) =  t[2];
	__LH_M4E(m,2,2) = -f[2];
	__LH_M4E(m,2,3) =  0.f;

	__LH_M4E(m,3,0) =  0.f;
	__LH_M4E(m,3,1) =  0.f;
	__LH_M4E(m,3,2) =  0.f;
	__LH_M4E(m,3,3) =  1.f;

	mat4x4_translate_in_place(m, -eye[0], -eye[1], -eye[2]);
}

////////////////////////////////////////////////////////////////////////////////
// Quaternions type: quat / quat_t
////////////////////////////////////////////////////////////////////////////////

typedef vec4 quat;
typedef union {
	quat q;
	__LH_VECTOR_STRUCT(x,y,z,w);
} quat_t;

#ifndef LINMATH_H_NO_BUILDERS
	#define QUAT(x,y,z,w)     (quat){x,y,z,w}
	#define iQUAT_T(x,y,z,w)  {.v={x,y,z,w}}
	#define QUAT_T(x,y,z,w)   (quat_t)iQUAT_T(x,y,z,w)
	#define cQUAT_T(v4)       QUAT_T(v4[0],v4[1],v4[2],v4[3])
#endif

#define quat_norm vec4_norm

static inline void quat_identity(quat q)
{
	q[0] = q[1] = q[2] = 0.f;
	q[3] = 1.f;
}
static inline void quat_add(quat r, quat a, quat b)
{
	for(int i=0; i<4; ++i)
		r[i] = a[i] + b[i];
}
static inline void quat_sub(quat r, quat a, quat b)
{
	for(int i=0; i<4; ++i)
		r[i] = a[i] - b[i];
}
static inline void quat_mul(quat r, quat p, quat q)
{
	vec3 w;
	vec3_mul_cross(r, p, q);
	vec3_scale(w, p, q[3]);
	vec3_add(r, r, w);
	vec3_scale(w, q, p[3]);
	vec3_add(r, r, w);
	r[3] = p[3]*q[3] - vec3_mul_inner(p, q);
}
static inline void quat_scale(quat r, quat v, float s)
{
	for(int i=0; i<4; ++i)
		r[i] = v[i] * s;
}
static inline float quat_inner_product(quat a, quat b)
{
	float p = 0.f;
	for(int i=0; i<4; ++i)
		p += b[i]*a[i];
	return p;
}
static inline void quat_conj(quat r, quat q)
{
	for(int i=0; i<3; ++i)
		r[i] = -q[i];
	r[3] = q[3];
}
static inline void quat_rotate(quat r, float angle, vec3 axis) {
	vec3 v;
	vec3_scale(v, axis, sinf(angle/2));
	for(int i=0; i<3; ++i)
		r[i] = v[i];
	r[3] = cosf(angle/2);
}
static inline void quat_mul_vec3(vec3 r, quat q, vec3 v)
{
	/*
	  Method by Fabian 'ryg' Giesen (of Farbrausch):
	    t = 2 * cross(q.xyz, v)
	    v' = v + q.w * t + cross(q.xyz, t)
	*/

	vec3 t;
	vec3 q_xyz = {q[0], q[1], q[2]};
	vec3 u = {q[0], q[1], q[2]};

	vec3_mul_cross(t, q_xyz, v);
	vec3_scale(t, t, 2);

	vec3_mul_cross(u, q_xyz, t);
	vec3_scale(t, t, q[3]);

	vec3_add(r, v, t);
	vec3_add(r, r, u);
}
static inline void mat4x4_from_quat(mat4x4 M, quat q)
{
	float a = q[3];
	float b = q[0];
	float c = q[1];
	float d = q[2];
	float a2 = a*a;
	float b2 = b*b;
	float c2 = c*c;
	float d2 = d*d;

	__LH_M4E(M,0,0) = a2 + b2 - c2 - d2;
	__LH_M4E(M,0,1) = 2.f * (b*c + a*d);
	__LH_M4E(M,0,2) = 2.f * (b*d - a*c);
	__LH_M4E(M,0,3) = 0.f;

	__LH_M4E(M,1,0) = 2.f * (b*c - a*d);
	__LH_M4E(M,1,1) = a2 - b2 + c2 - d2;
	__LH_M4E(M,1,2) = 2.f * (c*d + a*b);
	__LH_M4E(M,1,3) = 0.f;

	__LH_M4E(M,2,0) = 2.f * (b*d + a*c);
	__LH_M4E(M,2,1) = 2.f * (c*d - a*b);
	__LH_M4E(M,2,2) = a2 - b2 - c2 + d2;
	__LH_M4E(M,2,3) = 0.f;

	__LH_M4E(M,3,0) = __LH_M4E(M,3,1) = __LH_M4E(M,3,2) = 0.f;
	__LH_M4E(M,3,3) = 1.f;
}

static inline void mat4x4o_mul_quat(mat4x4 R, mat4x4 M, quat q)
{
	/*
	  NOTE: The way this is written only works for orthogonal matrices.
	  TODO: Take care of non-orthogonal case.
	*/

	vec4 m0, m1, m2;
	mat4x4_col(m0, M, 0);
	mat4x4_col(m1, M, 1);
	mat4x4_col(m2, M, 2);

	vec4 r0, r1, r2;
	quat_mul_vec3(r0, q, m0);
	quat_mul_vec3(r1, q, m1);
	quat_mul_vec3(r2, q, m2);

	mat4x4_set_col(R, r0, 0);
	mat4x4_set_col(R, r1, 1);
	mat4x4_set_col(R, r2, 2);
	__LH_M4E(R,3,0) = __LH_M4E(R,3,1) = __LH_M4E(R,3,2) = 0.f;
	__LH_M4E(R,3,3) = 1.f;
}
static inline void quat_from_mat4x4(quat q, mat4x4 M)
{
	float r = 0.f;
	int perm[] = {0, 1, 2, 0, 1};
	int *p = perm;

	for(int i=0; i<3; i++) {
		float m = M[i][i];
		if(m < r)
			continue;
		m = r;
		p = &perm[i];
	}

	r = sqrtf( 1.f + M[p[0]][p[0]] - M[p[1]][p[1]] - M[p[2]][p[2]] );

	if(r < LINMATH_H_NEAR_ZERO) {
		q[0] = 1.f;
		q[1] = q[2] = q[3] = 0.f;
		return;
	}

	q[0] = r/2.f;
	q[1] = ( __LH_M4E(M,p[0],p[1]) - __LH_M4E(M,p[1],p[0]) ) / (2.f*r);
	q[2] = ( __LH_M4E(M,p[2],p[0]) - __LH_M4E(M,p[0],p[2]) ) / (2.f*r);
	q[3] = ( __LH_M4E(M,p[2],p[1]) - __LH_M4E(M,p[1],p[2]) ) / (2.f*r);
}

#undef __LH_M4E
#undef __LH_POS
#undef __LH_NEG

#endif // LINMATH_H
