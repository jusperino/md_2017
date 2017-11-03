#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .d, .o and executable files
#

# define the executable file 
MAIN = blatt1

# define the C source files
SRCS = blatt1.cpp world.cpp timediscretization.cpp gravitypotential.cpp velocityverlet.cpp observer.cpp

# define the C compiler to use
CC = g++

# define any compile-time flags
CFLAGS = -Wall -g

# define any directories containing header files other than /usr/include
#
INCLUDES =

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS =

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link libm.so):
LIBS = -lm

SUFFIX=.cpp

# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:$(SUFFIX)=.o)

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above 
#

.PHONY: depend clean

all:    $(MAIN)

doc:	$(SRCS)
	doxygen

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
$(SUFFIX).o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@
	$(CC) -MM $(CFLAGS) $(INCLUDES) $< > $*.d

clean:
	$(RM) *.o *~ $(MAIN) *.d
	$(RM) -r api-doc

