# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 

# Source and object files
SRCS = main.cpp EdgeDisplacement.cpp ScrewDisplacement.cpp ParseLAMMPS.cpp ParseVASP.cpp
OBJS = $(SRCS:.cpp=.o)

# Output executable name
TARGET = run

# Default target
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Rule to compile .cpp to .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJS) $(TARGET)
