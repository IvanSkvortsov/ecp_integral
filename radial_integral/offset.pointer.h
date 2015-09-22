#ifndef __OFFSET_POINTER_H__
#define __OFFSET_POINTER_H__

// offset_pointer
//
// contains two pointers : _data, _offset
// the aim that I want to reach is fast usage of arrays which don't start with zero element array[0], 
//   e.i. from array[-4] to array[5] (size = 10) instead of array[0] to array[9]
template<class T>
class offset_pointer
{
	typedef T            value_type;
	typedef T *             pointer;
	typedef T const * const_pointer;
protected:
	T * _data, * _offset;
	void set_zero()
	{
		this->_data = 0;
		this->_offset = 0;
	}
public:
	offset_pointer():_data(0), _offset(0){}
	offset_pointer(T * __data):_data( __data ), _offset(_data){}
	offset_pointer(T * __data, T * __offset):_data( __data ), _offset(__offset){}
	offset_pointer(T * __data, int __offset):_data( __data ), _offset(_data + __offset){}
	~offset_pointer(){this->set_zero();}
	offset_pointer<T> & operator=(T * _ptr)
	{
		_data = _ptr;
		_offset = _data;
		return *this;
	}
	// data()
	void data(pointer __data){_data = __data;}
	pointer & data(){return _data;}
	const_pointer & data()const{return _data;}
	// offset()
	void offset(int __offset){_offset = (_data ? _data + __offset : _data);}
	void offset(pointer __offset){_offset = __offset;}
	pointer & offset(){return _offset;}
	const_pointer & offset()const{return _offset;}
	// operator[]
	T & operator[](int i){return _data[i];}
	T const & operator[](int i)const{return _data[i];}
	// operator()
	T & operator()(int i){return _offset[i];}
	T const & operator()(int i)const{return _offset[i];}
};

#endif//__OFFSET_POINTER_H__
