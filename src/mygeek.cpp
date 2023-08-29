#include <iostream>
#include "mygeek.h"

MyGeek::MyGeek(int val) : value(val) {}

void MyGeek::HelloGeek() {
    std::cout << "Hello geek! The value is: " << value << std::endl;
}