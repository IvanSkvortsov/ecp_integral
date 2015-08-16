#ifndef __OMEGA_INTEGRAL_H__
#define __OMEGA_INTEGRAL_H__
#include<cmath>// abs
#include"angular.omega.xyz.h"
#include"../spherical_harmonics/spherical.h"
#include<vector>

template<class T>
T omega_sph(spherical<T> const & sph, angular_omega_xyz<T> const & omega_xyz_)
{
	T value = T(0);
	typename spherical<T>::polynomial_type const * p = &sph[0];
	for(int i = 0; i < sph.size(); ++i)
	{
		if( !p->is_even_x() ) continue;
		value += p->d * omega_xyz_(p->x, p->y, p->z);
	}
	return value;
}

template<class T>
T omega_xyz_sph(int x, int y, int z, spherical<T> const & sph, angular_omega_xyz<T> const & omega_xyz_)
{
	typename spherical<T>::polynomial_type xyz(T(1), x, y, z), pol = xyz;
	typename spherical<T>::polynomial_type const * p = &sph[0];
	T value = T(0);
	for(int i = 0; i < sph.size(); ++i)
	{
		pol = *p * xyz;
		if( !pol.is_even_x() ) continue;
		value += pol.d * omega_xyz_(pol.x, pol.y, pol.z);
	}
	return value;
}

template<class T>
T omega_xyz_sph2(int x, int y, int z, spherical<T> const & sph1, spherical<T> const & sph2, angular_omega_xyz<T> const & omega_xyz_,
		spherical<T> & sph_buf)
{
	sph_buf  = sph1;
	sph_buf *= sph2;
	sph_buf.optimize_ez();
	sph_buf *= (typename spherical<T>::polynomial_type(T(1), x, y, z));
	return omega_sph<T>( sph_buf, omega_xyz_ );
}

bool check_xyz(int const & x, int const & y, int const & z)
{
	if( x < 0 || y < 0 || z < 0 ) return 1;
	return 0;
}

bool check_lm(int const & l, int const & m)
{
	if( l < 0 || abs(m) > l ) return 1;
	return 0;
}

template<class T>
T omega_xyz_sph2(int x, int y, int z, int l_a, int m_a, int l_b, int m_b, std::vector<std::vector<spherical<T> > > const & v2sph,
		angular_omega_xyz<T> const & omega_xyz_, spherical<T> & sph_buf)
{
	if( check_xyz( x, y, z ) || check_lm(l_a, m_a) || check_lm(l_b, m_b) )
	{
		std::cerr << "Error: [omega_xyz_sph2<T>(int x, int y, int z, int la, int ma, int lb, int mb, ...)]" << std::endl;
		std::cerr << "x  : " << x << std::endl;
		std::cerr << "y  : " << y << std::endl;
		std::cerr << "z  : " << z << std::endl;
		std::cerr << "la : " << l_a << std::endl;
		std::cerr << "ma : " << m_a << std::endl;
		std::cerr << "lb : " << l_b << std::endl;
		std::cerr << "mb : " << m_b << std::endl;
		exit(1);
	}
	int sum = x + y + z, max_l = l_a > l_b ? l_a : l_b, min_l = l_a > l_b ? l_b : l_a;
	if( sum + min_l < max_l ) return T(0);
	if( !sum && (l_a != l_b || m_a != m_b) ) return T(0);
	return omega_xyz_sph2<T>(x, y, z, v2sph[l_a][l_a + m_a], v2sph[l_b][l_b + m_b], omega_xyz_, sph_buf);
}

//-----------------------------------------------------------------------------------------------------------//

template<class T>
T omega_xyz_sph_l(int x, int y, int z, int l_a, int m_a, int l_b, T const * r, std::vector<std::vector<spherical<T> > > const & v2sph,
		angular_omega_xyz<T> const & omega_xyz_, spherical<T> & sph_buf)
{
	if( check_xyz( x, y, z ) || check_lm(l_a, m_a) || l_b < 0 || !r )
	{
		std::cerr << "Error: [omega_xyz_sph_l<T>(int x, int y, int z, int la, int lb, int lc, T r[3], ...)]" << std::endl;
		std::cerr << "x  : " << x << std::endl;
		std::cerr << "y  : " << y << std::endl;
		std::cerr << "z  : " << z << std::endl;
		std::cerr << "la : " << l_a << std::endl;
		std::cerr << "ma : " << m_a << std::endl;
		std::cerr << "lb : " << l_b << std::endl;
		std::cerr << "r  : [" << r << "]" << std::endl;
		exit(1);
	}
	int sum = x + y + z;
	if( !sum && l_a != l_b ) return T(0);
	T value = T(0);
	spherical<T> const & sph_a = v2sph[l_a][l_a + m_a];
	std::vector<spherical<T> > const & vsph_b = v2sph[l_b];
	for(int mb = -l_b; mb <= l_b; ++mb)
	{
		value += omega_xyz_sph2<T>(x, y, z, sph_a, vsph_b[l_b + mb], omega_xyz_, sph_buf) * vsph_b[l_b + mb].calc_r( r );
	}
	return value;
}

//-----------------------------------------------------------------------------------------------------------//

struct omega_index
{
	int a, b, c, d, e, f, l, lmb1, lmb2;
	omega_index():a(0), b(0), c(0), d(0), e(0), f(0), l(0), lmb1(0), lmb2(0){}
	omega_index(int const * A, int const * B):a(*A), b(*(A+1)), c(*(A+2)), d(*B), e(*(B+1)), f(*(B+2)), l(0), lmb1(0), lmb2(0){}
	//omega_index(int const & __a, int const & __b, int const & __c):a(__a), b(__b), c(__c), d(0), e(0), f(0), l(0), lmb1(0), lmb2(0){}
	void set_indexA(int const & __a, int const & __b, int const & __c)
	{
		a = __a; b = __b; c = __c;
	}
	void set_indexB(int const & __d, int const & __e, int const & __f)
	{
		d = __d; e = __e; f = __f;
	}
	int sum_A()const{return a + b + c;}
	int sum_B()const{return d + e + f;}
	int max_l()const{return (l > lmb1 ? (l > lmb2 ? l : lmb2) : (lmb1 > lmb2 ? lmb1 : lmb2));}
};

template<class T>
T omega_sph3(omega_index const & idx, T const * r_1, T const * r_2, std::vector<std::vector<spherical<T> > > const & v2sph,
		angular_omega_xyz<T> const & omega_xyz_, spherical<T> & sph_buf)
{
	if( check_xyz( idx.a, idx.b, idx.c ) || check_xyz( idx.d, idx.e, idx.f ) || idx.l < 0 || idx.lmb1 < 0 || idx.lmb2 < 0 || !r_1 || !r_2 )
	{
		std::cerr << "Error: [omega_xyz_sph3<T>(int x, int y, int z, int la, int lb, int lc, T rb[3], T rc[3], ...)]" << std::endl;
		std::cerr << "a    : " << idx.a << std::endl;
		std::cerr << "b    : " << idx.b << std::endl;
		std::cerr << "c    : " << idx.c << std::endl;
		std::cerr << "d    : " << idx.d << std::endl;
		std::cerr << "e    : " << idx.e << std::endl;
		std::cerr << "f    : " << idx.f << std::endl;
		std::cerr << "l    : " << idx.l << std::endl;
		std::cerr << "lmb1 : " << idx.lmb1 << std::endl;
		std::cerr << "lmb2 : " << idx.lmb2 << std::endl;
		std::cerr << "r1 : [" << r_1 << "]" << std::endl;
		std::cerr << "r2 : [" << r_2 << "]" << std::endl;
		exit(1);
	}
	int sum_a = idx.sum_A(), sum_b = idx.sum_B();
	if( sum_a > idx.lmb1 || sum_b > idx.lmb2 )
	{
		std::cerr << "Error: [omega_xyz_sph3<T>(int x, int y, int z, int la, int lb, int lc, T rb[3], T rc[3], ...)]" << std::endl;
		std::cerr << "sum_a : " << sum_a << std::endl;
		std::cerr << "lmb_a : " << idx.lmb1 << std::endl;
		std::cerr << "sum_b : " << sum_b << std::endl;
		std::cerr << "lmb_b : " << idx.lmb2 << std::endl;
		exit(1);
	}
	if( !sum_a && idx.l != idx.lmb1 ) return T(0);
	if( !sum_b && idx.l != idx.lmb2 ) return T(0);
	if( (sum_a + sum_b + idx.lmb1 + idx.lmb2)%2 == 0 && abs(idx.lmb1-idx.lmb2) > sum_a + sum_b )
		return T(0);
	T value = T(0), _omega_a = T(0), _omega_b = T(0);
	for(int m = -idx.l; m <= idx.l; ++m)
	{
		_omega_a = omega_xyz_sph_l<T>(idx.a, idx.b, idx.c, idx.l, m, idx.lmb1, r_1, v2sph, omega_xyz_, sph_buf);
		_omega_b = omega_xyz_sph_l<T>(idx.d, idx.e, idx.f, idx.l, m, idx.lmb2, r_2, v2sph, omega_xyz_, sph_buf);
		value += _omega_a * _omega_b;
	}
	return value;
}

#endif//__OMEGA_INTEGRAL_H__
