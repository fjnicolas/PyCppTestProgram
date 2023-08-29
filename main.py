import mygeek

geek = pybind_wrapper.MyGeek(2)

print("Result from C++ function:", geek.HelloGeek())
