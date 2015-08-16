#ifndef __ANGULAR_OMEGA_XYZ_H__
#define __ANGULAR_OMEGA_XYZ_H__
#include<cstddef>// size_t
#include<cstdlib>// exit
//#include"../../../lib_math/math.functions.h"// factt
#include"../../lib_math/math.functions.h"// factt

#ifndef _100x1024x1024_// 100 Mb
#define _100x1024x1024_ 104857600
#endif

#ifndef _1024x1024_// 1 Mb
#define _1024x1024_ 1048576
#endif

#ifndef _512x1024_// 512 Kb
#define _512x1024_  524288
#endif

#ifndef _256x1024_// 256 Kb
#define _256x1024_  262144
#endif

#ifndef _128x1024_// 128 Kb
#define _128x1024_  131072
#endif

// error
#define __ERROR_angular_omega_xyz__
#ifdef  __ERROR_angular_omega_xyz__
#include<iostream>
#endif


// info
#define __INFO_angular_omega_xyz__
#ifdef  __INFO_angular_omega_xyz__
#include<iostream>
#include<iomanip>
#endif

// login
#define __LOG_angular_omega_xyz__
#ifdef  __LOG_angular_omega_xyz__
#include<iostream>
#endif

template<class T, std::size_t _size = std::size_t(_1024x1024_)>
struct angular_omega_xyz
{
	typedef T value_type;
	typedef angular_omega_xyz<T, _size> omega_type;
protected:
	std::size_t _i_max, _j_max, _k_max;
	std::size_t _half_imax__p1, _half_jmax__p1, _half_kmax__p1;
	char _data[_size];
	T * _begin, * _end;
public:
	//
	angular_omega_xyz():
		_i_max(0), _j_max(0), _k_max(0), 
		_begin((T *)_data), _end( (T *)_data + _size / sizeof(T) ), 
		_half_imax__p1(0), _half_jmax__p1(0), _half_kmax__p1(0)
	{
#ifdef  __LOG_angular_omega_xyz__
		log("angular_omega_xyz()");
#endif
		T * _ptr = _begin;
		for(int i = 0; i < this->max_size(); ++i)
			new (_ptr++) T();
	}
	~angular_omega_xyz()
	{
#ifdef  __LOG_angular_omega_xyz__
		log("~angular_omega_xyz()");
#endif
	}
	std::size_t const & i_max()const{return _i_max;}
	std::size_t const & j_max()const{return _j_max;}
	std::size_t const & k_max()const{return _k_max;}
	std::size_t & i_max(){return _i_max;}
	std::size_t & j_max(){return _j_max;}
	std::size_t & k_max(){return _k_max;}
	std::size_t i_max(std::size_t const & __i_max){_half_imax__p1 = __i_max/2 + 1;return _i_max = __i_max;}
	std::size_t j_max(std::size_t const & __j_max){_half_jmax__p1 = __j_max/2 + 1;return _j_max = __j_max;}
	std::size_t k_max(std::size_t const & __k_max){_half_kmax__p1 = __k_max/2 + 1;return _k_max = __k_max;}
	std::size_t const max_size()const{return _end-_begin;}
	std::size_t const size()const{return _half_imax__p1 * _half_jmax__p1 * _half_kmax__p1;} 
	T * begin(){return _begin;}
	T const * begin()const{return _begin;}
	T * end(){return _end;}
	T const * end()const{return _end;}
	void set_zero()
	{
		T * p = _begin;
		for(std::size_t i = 0; i < this->size(); ++i)
			*p++ = T(0);
	}
	void run(std::size_t const & __i_max, std::size_t const & __j_max, std::size_t const & __k_max)
	{
		this->i_max(__i_max);
		this->j_max(__j_max);
		this->k_max(__k_max);
		this->run();
	}
	void run()
	{
		if( this->size() > this->max_size() )
		{
#ifdef  __ERROR_angular_omega_xyz__
			std::cerr << "Error: angular_omega_xyz<T>::run() : maximum size reached" << std::endl;
			std::cerr << "size     : " << this->size() << std::endl;
			std::cerr << "max_size : " << this->max_size() << std::endl;
#endif
			exit(1);
		}
		T * p = _begin;
		// TODO: 1) realize usage of _data empty space in run() method (to avoid memory allocation for pointers)
		// TODO: 2) realize less calculation and less mrmory usage for integral storage because of equivalence of (i,j,k), (j,i,k), etc.
		std::size_t __ijk_max = _i_max + _j_max + _k_max, __i_p_j;
		T * f_im1 = new T[_i_max/2+1], * f_jm1 = new T[_half_jmax__p1], * f_km1 = new T[_half_kmax__p1], * f_ijk_p1 = new T[__ijk_max/2+1];
		//T * f_im1 = new T[_i_max/2+1], * f_jm1 = new T[_j_max/2+1], * f_km1 = new T[_k_max/2+1], * f_ijk_p1 = new T[__ijk_max/2+1];
		T * p_im1 = f_im1, * p_jm1 = f_jm1, * p_km1 = f_km1;
		T f_im1_x_jm1;
		run_f_m1(f_im1, 0, _i_max);// (i - 1)!!
		run_f_m1(f_jm1, 0, _j_max);// (j - 1)!!
		run_f_m1(f_km1, 0, _k_max);// (k - 1)!!
		run_f_m1(f_ijk_p1, 2, __ijk_max+2);// (i + j + k + 1)!!
		for(std::size_t i = 0; i <= _i_max; i+=2)
		{
			p_jm1 = f_jm1;
			for(std::size_t j = 0; j <= _j_max; j+=2)
			{
				f_im1_x_jm1 = *p_im1 * *p_jm1;
				__i_p_j = i + j;
				p_km1 = f_km1;
				for(std::size_t k = 0; k <= _k_max; k+=2)
				{
					//*p++ = f_ijk_p1[(__i_p_j + k + 1)/2] / ( f_im1_x_jm1 * *p_km1++ );
					//*p++ = ( f_im1_x_jm1 * *p_km1++) / f_ijk_p1[(__i_p_j + k + 1)/2];
					//?? __i_p_j even, k even, so their sum also even => int((sum + 1)/2) -eq int(sum/2)
					*p++ = ( f_im1_x_jm1 * *p_km1++) / f_ijk_p1[(__i_p_j + k)/2];
				}
				++p_jm1;
			}
			++p_im1;
		}
		delete [] f_im1;
		delete [] f_jm1;
		delete [] f_km1;
		delete [] f_ijk_p1;
	}
	T operator()(std::size_t i, std::size_t j, std::size_t k)const
	{
		if( i%2 == 1 || j%2 == 1 || k%2 == 1 ) return T(0);
		if( i > _i_max || j > _j_max || k > _k_max )
		{
#ifdef  __ERROR_angular_omega_xyz__
			std::cerr << "Error: operator()(std::size_t i, std::size_t j, std::size_t k)" << std::endl;
			std::cerr << "i_max    : " << this->i_max() << std::endl;
			std::cerr << "j_max    : " << this->j_max() << std::endl;
			std::cerr << "k_max    : " << this->k_max() << std::endl;
			std::cerr << "i        : " << i << std::endl;
			std::cerr << "j        : " << j << std::endl;
			std::cerr << "k        : " << k << std::endl;
			std::cerr << "size     : " << this->size() << std::endl;
			std::cerr << "max_size : " << this->max_size() << std::endl;
#endif
			exit(1);
		}
		return _begin[ ((i * _half_jmax__p1 + j) * _half_kmax__p1 + k)/2 ];
		//return _begin[ ((i * _half_jmax__p1 + j) * _half_kmax__p1 + k)>>1 ];
		//return _begin[ ((i * (_j_max/2+1)   + j) * (_k_max/2+1)   + k)/2 ];
	}
	void pub_info()const
	{
#ifdef  __INFO_angular_omega_xyz__
		this->info();
#endif
	}
private:
	angular_omega_xyz(angular_omega_xyz<T> & );
	angular_omega_xyz(angular_omega_xyz<T> const & );
	angular_omega_xyz<T> & operator=(angular_omega_xyz<T> & );
	angular_omega_xyz<T> & operator=(angular_omega_xyz<T> const & );
	void run_f_m1(T * p, std::size_t __z_min, std::size_t __z_max)
	{
		for(std::size_t i = __z_min; i <= __z_max; i+=2)
			*p++ = factt<T>(i - 1);
	}
#ifdef  __INFO_angular_omega_xyz__
	void info()const
	{
		std::string s = "-----";
		int s_size = s.size();
		char * ptr_zero = 0x0;
		int addr = (char *)this - ptr_zero;
		int _pow = logi(addr, 16), leng = (s_size + 2) * 2 + (_pow + 3);
		int idx_w = 12, w = leng - idx_w;
		std::cout << s << " [" << this << "] " << s << std::endl;
		std::cout << std::setw(idx_w) << "max_size : " << std::setw(w) << this->max_size() << std::endl;
		std::cout << std::setw(idx_w) << "size : " << std::setw(w) << this->size() << std::endl;
		std::cout << std::setw(idx_w) << "i_max : " << std::setw(w) << this->_i_max << std::endl;
		std::cout << std::setw(idx_w) << "j_max : " << std::setw(w) << this->_j_max << std::endl;
		std::cout << std::setw(idx_w) << "k_max : " << std::setw(w) << this->_k_max << std::endl;
		for(int i = 0; i < leng; ++i) std::cout << '-'; std::cout << std::endl;
	}
#endif
#ifdef  __LOG_angular_omega_xyz__
	void log(std::string const & s)const
	{
		std::cout << "[" << this << "] angular_omega_xyz<T>::" << s << std::endl;
	}
#endif
};

#endif//__ANGULAR_OMEGA_XYZ_H__