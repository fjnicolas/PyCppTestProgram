from pybind_wrapper import MyGeek

geek = MyGeek(2)

print("Result from C++ function:", geek.HelloGeek())
