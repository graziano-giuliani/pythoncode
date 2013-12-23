
F2PY3 = f2py
F2PY2 = f2py2

F2PY3 = f2py --f90flags="-O0 -g -fbacktrace \
  -fdump-core -Wall -pedantic -Wconversion -Wsurprising -Waliasing \
  -Wampersand -Warray-bounds -Wcharacter-truncation \
  -Wline-truncation -Wintrinsics-std -Wno-tabs -Wunderflow \
  -Wunused-parameter -Wno-align-commons -fbounds-check -I."
F2PY2 = f2py2 --f90flags="-O0 -g -fbacktrace \
  -fdump-core -Wall -pedantic -Wconversion -Wsurprising -Waliasing \
  -Wampersand -Warray-bounds -Wcharacter-truncation \
  -Wline-truncation -Wintrinsics-std -Wno-tabs -Wunderflow \
  -Wunused-parameter -Wno-align-commons -fbounds-check -I."

all:: ossir.cpython-33m.so ossir.so

ossir.cpython-33m.so: oss_ir.F90
	$(F2PY3) -c -m ossir $<

ossir.so: oss_ir.F90
	$(F2PY2) -c -m ossir $<

clean:
	rm -f ossir*.so
