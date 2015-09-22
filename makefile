
ecp.integral.exe: ecp.integral.h ecp.integral.demo.h ecp.integral.main.cpp
	g++ ecp.integral.main.cpp radial_integral/recursion/spec_func/erfh.o radial_integral/recursion/spec_func/dawson.o -o $@
	./$@ ecp.integral.inp

clean:
	rm ecp.integral.exe
