UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
	CC=c++
	CFLAGS=-O3
	CPPFLAGS=-std=c++14 -I/usr/local/include
	LDFLAGS=-L/usr/local/lib 
	LDFLAGSGUI=-L/usr/local/lib -framework Cocoa -framework OpenGL -framework IOKit -framework CoreVideo -lGLEW -lglfw3
endif
ifeq ($(UNAME), Linux)
	CC=g++
	CFLAGS=-O3
	CPPFLAGS=-std=c++14 -I/usr/local/include
	LDFLAGS=-L/usr/local/lib
	LDFLAGSGUI=-L/usr/local/lib  -lGL -lglfw -lGLEW
	# -lX11 -lXxf86vm -lXrandr -lpthread -lXi -ldl -lXinerama -lXcursor
endif

%: %.cc
	$(CC) $< -o $@ $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) 

demo2gui: demo2gui.cc
	$(CC) $< -o $@ $(CFLAGS) $(CPPFLAGS) $(LDFLAGSGUI) 
	
all: demo1text demo2gui main

clean:
	rm -Rf DATA
	rm demo1text demo2gui main
	
