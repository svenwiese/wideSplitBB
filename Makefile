###################
# modify paths... #
###################

#on linux...
#BASE=/opt/ibm/ILOG/CPLEX_Studio126/cplex
#DEVCPLEXLIBDIR	= ${BASE}/lib/x86-64_linux/static_pic

#on mac...
BASE=/Users/sven/Applications/IBM/ILOG/CPLEX_Studio1261/cplex
DEVCPLEXLIBDIR	= ${BASE}/lib/x86-64_osx/static_pic

#---------------------------------
DEVCPLEXINCDIR	= ${BASE}/include/ilcplex/

CXX		= g++
FLAGS		=
CFLAGS		= -I$(DEVCPLEXINCDIR) -O3 -DNDEBUG -pipe -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings
LDFLAGS		= -L$(DEVCPLEXLIBDIR) -lcplex -lpthread -lm

SRCDIR		= src
OBJDIR		= obj
TARGET		= branch_and_hole
OBJ		= branch_and_hole.o holes.o

SRCFILES	= $(addprefix $(SRCDIR)/,$(OBJ:.o=.cpp))
OBJFILES	= $(addprefix $(OBJDIR)/,$(OBJ))

$(TARGET): $(OBJFILES)
		$(CXX) $(FLAGS) $(OBJFILES) $(LDFLAGS) -o $@

.PHONY:	clean
clean:
		rm -f $(OBJFILES) $(TARGET)

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp $(OBJDIR)
		$(CXX) $(FLAGS) $(CFLAGS) -c $< -o $@

all: ${OBJFILES} ${TARGET}

.PHONY: zip
zip:
	mkdir -p lib
	mkdir -p include
	mkdir -p include/ilcplex
	cp $(DEVCPLEXLIBDIR)/lib*.a lib/
	cp $(DEVCPLEXINCDIR)/ilcplex/cplex.h include/ilcplex/
	zip watsoncplex.zip Makefile src/* lib/* include/ilcplex/*
