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

//---------------------------------------------------------------------------------------//
//----------------------------------           ------------------------------------------//
//---------------------------------  OMEGA_INT  -----------------------------------------//
//----------------------------------           ------------------------------------------//
//---------------------------------------------------------------------------------------//
template<class T>
T omega_int(typename spherical<T>::polynomial_type const & pol, spherical<T> const & sph_x, spherical<T> const & sph_pp, 
		angular_omega_xyz<T> const & omega_xyz_, spherical<T> & sph_buf)
{
	int sum_x = pol.sum_x(), n_x = sph_x.n(), n_pp = sph_pp.n(), m_x = sph_x.m(), m_pp = sph_pp.m();
	int n_max = n_x > n_pp ? n_x : n_pp, n_min = n_x + n_pp - n_max;
	if( !sum_x && (n_x != n_pp || m_x != m_pp) ) return T(0);
	if( sum_x + n_min < n_max ) return T(0);
	sph_buf  = sph_x;
	sph_buf *= sph_pp;
	sph_buf.optimize_ez();
	sph_buf *= pol;
	T value = T(0);
	typename spherical<T>::polynomial_type const * p = &sph_buf[0];
	for(int i = 0; i < sph_buf.size(); ++i)
	{
		if( !p->is_even_x() ) continue;
		value += p->d * omega_xyz_(p->x, p->y, p->z);
	}
	return value;
}

template<class T>
T omega_int(typename spherical<T>::polynomial_type const & pol, std::vector<T> const & sph_r, std::vector<spherical<T> > const & vsph_x,
		spherical<T> const & sph_pp, angular_omega_xyz<T> const & omega_xyz_, spherical<T> & sph_buf)
{
	T value = T(0);
	if( sph_r.size() != vsph_x.size() )
	{
		std::cerr << "Error: omega_int(polynom pol, vector sph_r, vector<spherical> vsph, spherical sph_pp, omega_xyz, sph_buf)" << std::endl;
		std::cerr << "sph_r.size : " << sph_r.size() << std::endl;
		std::cerr << " vsph.size : " << vsph_x.size() << std::endl;
		exit(1);
	}
	for(int i = 0; i < vsph_x.size(); ++i)
		value += sph_r[i] * omega_int<T>( pol, vsph_x[i], sph_pp, omega_xyz_, sph_buf);
	return value;
}

int check_xyz(int const * x)
{
	if( x[0] < 0 || x[1] < 0 || x[2] < 0 )
		return 1;
	return 0;
}

template<class T>
T omega_int(    typename spherical<T>::polynomial_type const & pol_a, typename spherical<T>::polynomial_type const & pol_b,
		std::vector<T> const & sph_ra, std::vector<T> const & sph_rb,
		std::vector<spherical<T> > const & vsph_a, std::vector<spherical<T> > const & vsph_b,
		std::vector<spherical<T> > const & vsph_pp,
		angular_omega_xyz<T> const & omega_xyz_, spherical<T> & sph_buf)
{
	if( check_xyz( &pol_a.x ) || check_xyz( &pol_b.x ) )
	{
		std::cerr << "Error: [omega_int(polynomial const & pol_a, polynomial const & pol_b, ...)]" << std::endl;
		std::cerr << "pol_a : " << pol_a.x << ' ' << pol_a.y << ' ' << pol_a.z << std::endl;
		std::cerr << "pol_b : " << pol_b.x << ' ' << pol_b.y << ' ' << pol_b.z << std::endl;
		exit(1);
	}
	int l_pp = vsph_pp[0].n();
	int sum_b = pol_b.sum_x(), l_b = vsph_b[0].n();
	int l_max = l_b > l_pp ? l_b : l_pp, l_min = l_b + l_pp - l_max;
	if( !sum_b && l_b != l_pp ) return T(0);
	if( sum_b + l_min < l_max ) return T(0);
	//
	int sum_a = pol_a.sum_x(), l_a = vsph_a[0].n();
	l_max = l_a > l_pp ? l_a : l_pp;
	l_min = l_a + l_pp - l_max;
	if( !sum_a && l_a != l_pp ) return T(0);
	if( sum_a + l_min < l_max ) return T(0);
	//
	if( (sum_a + sum_b + l_a + l_b)%2==0 && abs(l_a - l_b) > (sum_a + sum_b) ) return T(0);
	//
	T value = T(0);
	for(int i = 0; i < vsph_pp.size(); ++i)
		value += omega_int<T>( pol_a, sph_ra, vsph_a, vsph_pp[i], omega_xyz_, sph_buf ) *
			 omega_int<T>( pol_b, sph_rb, vsph_b, vsph_pp[i], omega_xyz_, sph_buf );
	return value;
}

//-------------------------------------------//
//-----------                  --------------//
//---------- Newton coefficient -------------//
//-----------                  --------------//
//-------------------------------------------//
template<class T, class vector_type>
void run_nc(vector_type & v, std::size_t n)
{
	v.resize(n + 1);
	int i = 0;
	for(i = 0; i < v.size()/2; ++i)
	{
		v[i] = NewtonC<T>(n, i);
		v[n-i] = v[i];
	}
	if(v.size()%2) v[i] = NewtonC<T>(n, i);
}

template<class T, class vector_type>
void run_cx(vector_type & v, T const & R, std::size_t n)
{
	v.resize(n + 1);
	v[0] = T(1);
	for(int i = 1; i < v.size(); ++i)
		v[i] = v[i-1] * R;
}

template<class T>
struct newtc_x_ca
{
	std::vector<T> nca_x, nca_y, nca_z;
	void run(T const * R, std::size_t const * ax)
	{
		run(R, ax[0], ax[1], ax[2]);
	}
	void run(T const * R, std::size_t ax, std::size_t ay, std::size_t az)
	{
		run_v(nca_x, R[0], ax);
		run_v(nca_y, R[1], ay);
		run_v(nca_z, R[2], az);
	}
	int n_x()const{return nca_x.size()-1;}
	int n_y()const{return nca_y.size()-1;}
	int n_z()const{return nca_z.size()-1;}
	T operator()(int i, int j, int k)const{return nca_x[i] * nca_y[j] * nca_z[k];}
	T operator()(int const * i)const{return nca_x[*i] * nca_y[ i[1] ] * nca_z[ i[2] ];}
protected:
	void run_v(std::vector<T> & v, T const & R, std::size_t n)
	{
		v.resize(n + 1);
		std::vector<T> nc_v, ca_v;
		run_nc<T, std::vector<T> >(nc_v, n);
		run_cx<T, std::vector<T> >(ca_v, R, n);
		for(int i = 0; i < v.size(); ++i)
			v[i] = nc_v[i] * ca_v[n-i];
	}
};

struct _abc_
{
	int a, b, c;
	_abc_():a(0), b(0), c(0){}
	_abc_(int __a, int __b, int __c):a(__a), b(__b), c(__c){}
	void set(int __a, int __b, int __c)
	{
		a = __a; b = __b; c = __c;
	}
	int sum()const{return a + b + c;}
	int & x(){return a;}
	int & y(){return b;}
	int & z(){return c;}
	int const & x()const{return a;}
	int const & y()const{return b;}
	int const & z()const{return c;}
};

template<class T>
int run_n_xyz(std::vector<T> & v, int n, int ax, int ay, int az)
{
	v.reserve(n * n + 1);
	v.clear();
	int nmi = 0;
	for(int i = 0; i <= ax && i <= n; ++i)
	{
		nmi = n - i;
		for(int j = 0; j <= ay && j <= nmi; ++j)
			v.push_back( T(i, j, nmi-j) );
	}
}

struct v2_xyz
{
	int l, x, y, z;
	std::vector<std::vector<_abc_> > _v2_xyz;
	v2_xyz():l(0), x(0), y(0), z(0){}
	void run(int __l, int __x, int __y, int __z)
	{
		if( __l < 0 || (__x + __y + __z) != __l || __x < 0 || __y < 0 || __z < 0 )
		{
			std::cerr << "Error: [v2_xyz<T>::run(int l, int x, int y, int z)]" << std::endl;
			std::cerr << "l : " << __l << std::endl;
			std::cerr << "x : " << __x << std::endl;
			std::cerr << "y : " << __y << std::endl;
			std::cerr << "z : " << __z << std::endl;
			exit(1);
		}
		l = __l; x = __x; y = __y; z = __z;
		_v2_xyz.resize(l + 1);
		for(int i = 0; i < _v2_xyz.size(); ++i)
			run_n_xyz<_abc_>(_v2_xyz[i], i, x, y, z);
	}
	std::size_t size()const{return _v2_xyz.size();}
	std::vector<_abc_> & operator[](std::size_t i){return _v2_xyz[i];}
	std::vector<_abc_> const & operator[](std::size_t i)const{return _v2_xyz[i];}
};

//---------------------------------------------------------//
//-------------------              ------------------------//
//------------------ OMEGA INTEGRAL -----------------------//
//-------------------              ------------------------//
//---------------------------------------------------------//

template<class T>
struct omega_integral
{
	std::vector<T> v;
	int la_sz, lb_sz, lmb_asz, lmb_bsz;
	omega_integral():v(std::size_t(_1024x1024_)/sizeof(T)), la_sz(0), lb_sz(0), lmb_asz(0), lmb_bsz(0)
	{
		for(int i = 0; i < v.size(); ++i) v[i] = T(0);
	}
	//
	T & operator[](std::size_t i){return v[i];}
	T const & operator[](std::size_t i)const{return v[i];}
	std::size_t max_size()const{return v.size();}
	std::size_t size()const{return la_sz * lb_sz * lmb_asz * lmb_bsz;}
	//
	std::size_t la_size()const{return la_sz;}
	std::size_t lb_size()const{return lb_sz;}
	std::size_t lmb_asize()const{return lmb_asz;}
	std::size_t lmb_bsize()const{return lmb_bsz;}
	//
	void run(int, _abc_ const &, _abc_ const & , T const *, T const *, angular_omega_xyz<T> const &);
	T & operator()(int na, int nb, int lmb_a, int lmb_b){return v[((na * lb_sz + nb) * lmb_asz + lmb_a) * lmb_bsz + lmb_b];}
	T const & operator()(int na, int nb, int lmb_a, int lmb_b)const{return v[((na * lb_sz + nb) * lmb_asz + lmb_a) * lmb_bsz + lmb_b];}
};

template<class T>
void norm_v3(T * norm_v, T const * v)
{
	T sqr_v = 0;
	for(int i = 0; i < 3; ++i) sqr_v += v[i] * v[i];
	if( sqr_v == T(0) )
		for(int i = 0; i < 3; ++i) norm_v[i] = v[i];
	else
	{
		sqr_v = sqrt( sqr_v );
		for(int i = 0; i < 3; ++i) norm_v[i] = v[i]/sqr_v;
	}
}

template<class T>
struct omega_integral_index
{
	std::vector<_abc_> const * p_xyz_a, * p_xyz_b;
	newtc_x_ca<T> const * p_nc_x_ca, * p_nc_x_cb;
	std::vector<T> const * sph_ra, * sph_rb;
	std::vector<spherical<T> > const * p_vsph_l, * p_vsph_lmb_a, * p_vsph_lmb_b;
	//
	omega_integral_index(): p_xyz_a(0), p_xyz_b(0), p_nc_x_ca(0), p_nc_x_cb(0), sph_ra(0), sph_rb(0),
				p_vsph_l(0), p_vsph_lmb_a(0), p_vsph_lmb_b(0){}
	// set
	void set_xyz_a(std::vector<_abc_> const & vxyz_a){p_xyz_a = &vxyz_a;}
	void set_xyz_b(std::vector<_abc_> const & vxyz_b){p_xyz_b = &vxyz_b;}
	void set_newt(newtc_x_ca<T> const & nc_x_ca, newtc_x_ca<T> const & nc_x_cb)
	{
		p_nc_x_ca = &nc_x_ca;
		p_nc_x_cb = &nc_x_cb;
	}
	void set_sph_ra(std::vector<T> const & s_ra){sph_ra = &s_ra;}
	void set_sph_rb(std::vector<T> const & s_rb){sph_rb = &s_rb;}
	void set_sph_l(std::vector<spherical<T> > const & vsph_l){p_vsph_l = &vsph_l;}
	void set_sph_lmb_a(std::vector<spherical<T> > const & vsph_lmb_a){p_vsph_lmb_a = &vsph_lmb_a;}
	void set_sph_lmb_b(std::vector<spherical<T> > const & vsph_lmb_b){p_vsph_lmb_b = &vsph_lmb_b;}
};

template<class T>
T run_omega(omega_integral_index<T> const & indx, angular_omega_xyz<T> const & omega_xyz_, spherical<T> sph_buf)
{
	std::vector<_abc_> const & vxyz_a = *indx.p_xyz_a, & vxyz_b = *indx.p_xyz_b;
	newtc_x_ca<T> const & nc_x_ca = *indx.p_nc_x_ca, & nc_x_cb = *indx.p_nc_x_cb;
	std::vector<T> const & sph_ra = *indx.sph_ra, & sph_rb = *indx.sph_rb;
	std::vector<spherical<T> > const & vsph_l = *indx.p_vsph_l, & vsph_lmb_a = *indx.p_vsph_lmb_a, & vsph_lmb_b = *indx.p_vsph_lmb_b;
	//
	typename spherical<T>::polynomial_type pol_a(T(1), 0, 0, 0), pol_b(T(1), 0, 0, 0);
	_abc_ const * xyz_a, * xyz_b;
	T newt_x_ca = 0, newt_x_cb = 0, value = 0;
	int l = indx.p_vsph_l->operator[](0).n();
	for(int i = 0; i < vxyz_a.size(); ++i)
	{
		xyz_a = &vxyz_a[i];
		pol_a.set_x( &xyz_a->x() );
		newt_x_ca = nc_x_ca(&xyz_a->x());
		for(int j = 0; j < vxyz_b.size(); ++j)
		{
			xyz_b = &vxyz_b[j];
			pol_b.set_x( &xyz_b->x() );
			newt_x_cb = nc_x_cb(&xyz_b->x());
			value += newt_x_ca * newt_x_cb * 
				omega_int<T>( pol_a, pol_b, sph_ra, sph_rb, vsph_lmb_a, vsph_lmb_b, vsph_l, omega_xyz_, sph_buf );
		}
	}
	return value;
}

template<class T>
void spherical_rx(std::vector<std::vector<T> > & v2sph_rx, std::vector<std::vector<spherical<T> > > const & v2_sph, T const * norm_rx)
{
	for(int i = 0; i < v2sph_rx.size(); ++i)
	{
		v2sph_rx[i].resize(2 * i + 1);
		for(int j = 0; j < v2sph_rx[i].size(); ++j)
			v2sph_rx[i][j] = v2_sph[i][j].calc_r(norm_rx);
	}
}

int max(int const & a, int const & b)
{
	return a > b ? a : b;
}

template<class T>
void omega_integral<T>::run(int l, _abc_ const & la, _abc_ const & lb, T const * ra, T const * rb, angular_omega_xyz<T> const & omega_xyz_)
{
	la_sz = la.sum() + 1;
	lb_sz = lb.sum() + 1;
	lmb_asz = l + la_sz;
	lmb_bsz = l + lb_sz;
	// 3d vector
	T norm_ra[3], norm_rb[3];
	norm_v3<T>( norm_ra, ra );
	norm_v3<T>( norm_rb, rb );
	// index
	v2_xyz v2_xyz_a, v2_xyz_b;
	v2_xyz_a.run( la.sum(), la.x(), la.y(), la.z() );
	v2_xyz_b.run( lb.sum(), lb.x(), lb.y(), lb.z() );
	// newtonc_x_ca
	newtc_x_ca<T> nc_x_ca, nc_x_cb;
	nc_x_ca.run(ra, la.x(), la.y(), la.z());
	nc_x_cb.run(rb, lb.x(), lb.y(), lb.z());
	// spherical
	std::vector<std::vector<spherical<T> > > v2_sph;
	std::vector<spherical<T> > * p_vsph = 0;
	spherical<T> sph_buf;
	sph_buf.reserve( 100 );
	std::cout << "lmb_max : " << l + max(la.sum(), lb.sum()) << std::endl;
	v2_sph.resize( l + max(la.sum(), lb.sum()) + 1);
	for(int i = 0; i < v2_sph.size(); ++i)
	{
		p_vsph = &v2_sph[i];
		p_vsph->resize(2 * i + 1);
		for(int j = 0; j < p_vsph->size(); ++j)
			(*p_vsph)[j].run(i, j-i);
	}
	// spherical r
	std::vector<std::vector<T> > v2sph_ra(l + la.sum() + 1), v2sph_rb(l + lb.sum() + 1);
	spherical_rx<T>( v2sph_ra, v2_sph, norm_ra );
	spherical_rx<T>( v2sph_rb, v2_sph, norm_rb );
	// omega_integral_index
	omega_integral_index<T> indx;
	indx.set_newt( nc_x_ca, nc_x_cb );
	indx.set_sph_l( v2_sph[l] );
	// run
	T * p = &v[0];
	for(int na = 0; na < la_sz; ++na)
	{
		indx.set_xyz_a( v2_xyz_a[na] );
		for(int nb = 0; nb < lb_sz; ++nb)
		{
			indx.set_xyz_b( v2_xyz_b[nb] );
			for(int lmb_a = 0; lmb_a < lmb_asz; ++lmb_a)
			{
				indx.set_sph_ra( v2sph_ra[lmb_a] );
				indx.set_sph_lmb_a( v2_sph[lmb_a] );
				for(int lmb_b = 0; lmb_b < lmb_bsz; ++lmb_b)
				{
					indx.set_sph_rb( v2sph_rb[lmb_b] );
					indx.set_sph_lmb_b( v2_sph[lmb_b] );
					*p++ = run_omega<T>(indx, omega_xyz_, sph_buf);
				}
			}
		}
	}
}


template<class T>
void print_omega_integral(std::ostream & out, omega_integral<T> const & o_i)
{
	int prec = 16, w = prec + 8;
	out.setf( std::ios::scientific );
	out.precision( prec );
	T const * p = &o_i[0];
	for(int na = 0; na < o_i.la_size(); ++na)
		for(int nb = 0; nb < o_i.lb_size(); ++nb)
			for(int lmb_a = 0; lmb_a < o_i.lmb_asize(); ++lmb_a)
				for(int lmb_b = 0; lmb_b < o_i.lmb_bsize(); ++lmb_b)
					out << std::setw(4) << na << std::setw(4) << nb << std::setw(4) << lmb_a << std::setw(4) << lmb_b <<
						std::setw( w ) << *p++ << std::endl;
}




#endif//__OMEGA_INTEGRAL_H__
