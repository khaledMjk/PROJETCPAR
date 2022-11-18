# include à modifier selon sa configuration
include make_msys.inc 


ALL= vortexSimulation.exe

default:	help
all: $(ALL)

clean:
	@rm -fr *.o *.exe *~ *.png

OBJS= vortex.o screen.o runge_kutta.o cloud_of_points.o cartesian_grid_of_speed.o vortexSimulation.o

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) -c -o $@ $^

vortexSimulation.exe: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIB)

help:
	@echo "Available targets : "
	@echo "    all                           : compile all executables"
	@echo "    syracuse_simple.exe           : compile simple syracuse sequence executable"
	@echo "    syracuse_orbite.exe           : compile syracuse with orbit executable"
	@echo "    simple_automate_cellulaire.exe: compile simple cellular automata executable"
	@echo "Add DEBUG=yes to compile in debug"
	@echo "Configuration :"
	@echo "    CXX      :    $(CXX)"
	@echo "    CXXFLAGS :    $(CXXFLAGS)"
